/*
 * grm_shard_tool -- bin + individual-block-jackknife a GRM that was computed
 * in row-chunked shards (e.g. `plink --make-grm-bin --parallel k n`),
 * without ever assembling the full dense .grm.bin.
 *
 * Subcommands:
 *   accumulate  -- read one shard, accumulate per-bin / per-jackknife-block
 *                  sums (not means) into a small output file
 *   merge       -- sum several shards' accumulator files and compute the
 *                  final bin means, SDs, and jackknife SEs
 *   ranges      -- print the row range a given shard *should* cover, for
 *                  manual sanity-checking (no shard file needed)
 *
 * The accumulation logic (full[]/drop[][] over Acc{sum,sum_sq,n}, and the
 * delete-block jackknife variance) is the same as full_grm_bin's
 * blocked_grm_bin_means.cpp, just split across the shard boundary: summing
 * Acc structs is associative, so accumulating per shard and adding the
 * results together afterward gives identical totals to one pass over the
 * whole matrix.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

struct ID {
    std::string fid;
    std::string iid;
    bool operator==(const ID& other) const {
        return fid == other.fid && iid == other.iid;
    }
};

struct IDHash {
    std::size_t operator()(const ID& x) const {
        return std::hash<std::string>()(x.fid) ^ (std::hash<std::string>()(x.iid) << 1);
    }
};

struct Bin {
    double left;
    double right;
};

struct Acc {
    double sum = 0.0;
    double sum_sq = 0.0;
    std::uint64_t n = 0;

    void add(double x) {
        sum += x;
        sum_sq += x * x;
        n += 1;
    }
    void add(const Acc& other) {
        sum += other.sum;
        sum_sq += other.sum_sq;
        n += other.n;
    }
};

// ---------------------------------------------------------------------------
// Shared readers (same formats as full_grm_bin/blocked_grm_bin_means.cpp)
// ---------------------------------------------------------------------------

static bool is_missing_token(const std::string& s) {
    return s == "NA" || s == "NaN" || s == "nan" || s == ".";
}

static double parse_double_or_nan(const std::string& s) {
    if (is_missing_token(s)) return std::numeric_limits<double>::quiet_NaN();
    try {
        std::size_t idx = 0;
        double x = std::stod(s, &idx);
        return idx == s.size() ? x : std::numeric_limits<double>::quiet_NaN();
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

static std::vector<ID> read_ids(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Could not open grm.id file: " + path);

    std::vector<ID> ids;
    std::string fid, iid;
    while (in >> fid >> iid) ids.push_back({fid, iid});
    if (ids.empty()) throw std::runtime_error("No IDs found in grm.id file: " + path);
    return ids;
}

static std::vector<double> read_pheno_aligned(const std::string& path, const std::vector<ID>& ids) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Could not open phenotype file: " + path);

    std::unordered_map<ID, double, IDHash> pheno;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string fid, iid, ystr;
        if (!(iss >> fid >> iid >> ystr)) continue;
        if (fid == "FID" && iid == "IID") continue;
        pheno[{fid, iid}] = parse_double_or_nan(ystr);
    }

    std::vector<double> y(ids.size(), std::numeric_limits<double>::quiet_NaN());
    for (std::size_t i = 0; i < ids.size(); ++i) {
        auto it = pheno.find(ids[i]);
        if (it != pheno.end()) y[i] = it->second;
    }
    return y;
}

static std::vector<Bin> read_bins(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Could not open bins file: " + path);

    std::vector<Bin> bins;
    double a, b;
    while (in >> a >> b) {
        if (!(a < b)) throw std::runtime_error("Invalid bin with left >= right");
        bins.push_back({a, b});
    }
    if (bins.empty()) throw std::runtime_error("No bins found in file: " + path);

    std::sort(bins.begin(), bins.end(), [](const Bin& x, const Bin& y) { return x.left < y.left; });
    return bins;
}

static int get_bin(double g, const std::vector<Bin>& bins) {
    int lo = 0, hi = static_cast<int>(bins.size()) - 1;
    while (lo <= hi) {
        int mid = lo + (hi - lo) / 2;
        const Bin& b = bins[mid];
        bool in_bin = (mid == static_cast<int>(bins.size()) - 1)
                          ? (g >= b.left && g <= b.right)
                          : (g >= b.left && g < b.right);
        if (in_bin) return mid;
        if (g < b.left) hi = mid - 1; else lo = mid + 1;
    }
    return -1;
}

// Deterministic given (n, nblocks, seed) alone -- every shard computes the
// same assignment independently, no coordination between shard jobs needed.
static std::vector<int> make_random_blocks(std::size_t n, int nblocks, unsigned int seed) {
    if (static_cast<std::size_t>(nblocks) > n) {
        throw std::runtime_error("Number of blocks exceeds number of individuals");
    }
    std::vector<std::size_t> perm(n);
    for (std::size_t i = 0; i < n; ++i) perm[i] = i;

    std::mt19937 rng(seed);
    std::shuffle(perm.begin(), perm.end(), rng);

    std::vector<int> block_of(n, -1);
    for (std::size_t rank = 0; rank < n; ++rank) {
        int b = static_cast<int>((static_cast<long double>(rank) * nblocks) / n);
        if (b >= nblocks) b = nblocks - 1;
        block_of[perm[rank]] = b;
    }
    return block_of;
}

static double bin_midpoint(const Bin& b) { return 0.5 * (b.left + b.right); }

static double safe_mean(double sum, std::uint64_t n) {
    return n == 0 ? std::numeric_limits<double>::quiet_NaN() : sum / static_cast<double>(n);
}

static double safe_sample_variance(double sum, double sum_sq, std::uint64_t n) {
    if (n < 2) return std::numeric_limits<double>::quiet_NaN();
    const double nd = static_cast<double>(n);
    double var = (sum_sq - (sum * sum) / nd) / static_cast<double>(n - 1);
    if (var < 0.0 && std::abs(var) < 1e-12) var = 0.0;
    return var;
}

static double safe_sample_sd(double sum, double sum_sq, std::uint64_t n) {
    double var = safe_sample_variance(sum, sum_sq, n);
    return std::isnan(var) ? var : std::sqrt(var);
}

static double safe_se_of_mean(double sum, double sum_sq, std::uint64_t n) {
    double var = safe_sample_variance(sum, sum_sq, n);
    return std::isnan(var) ? var : std::sqrt(var / static_cast<double>(n));
}

// ---------------------------------------------------------------------------
// Row-range recovery for one shard of `--parallel k n_shards`
// ---------------------------------------------------------------------------

// Cumulative number of GRM entries WRITTEN TO FILE (including the diagonal)
// in rows [0, i) -- confirmed against plink 1.9's calc_rel(): each row writes
// its off-diagonal entries followed by one diagonal (IBC) entry.
static uint64_t cumulative_entries(uint64_t i) { return i * (i + 1) / 2; }

struct ShardRange {
    uint64_t row_start;
    uint64_t row_end;  // exclusive
};

// This is plink 1.9's actual --parallel split point, not an approximation:
// see plink_common.c parallel_bounds()/triangle_divide() and calc_rel()'s
// `min_sample`/`max_parallel_sample` in plink_calc.c. Unlike the file layout
// above, the SPLIT BOUNDARY itself is chosen by balancing off-diagonal pair
// counts only (vv*(vv-1) >= target, not vv*(vv+1)) -- the diagonal is cheap
// to compute and just gets folded into whichever row it belongs to when the
// file is written. triangle_divide(0, -1) is defined to return 1; plink
// then snaps that 1 back to 0 so row 0 (which has no off-diagonal entries
// of its own) isn't dropped by the split.
static uint64_t triangle_divide_off_diag(uint64_t target) {
    if (target == 0) return 1;
    double approx = (std::sqrt(static_cast<double>(target)) + 1.0);
    int64_t v = static_cast<int64_t>(approx) + 1;
    if (v < 1) v = 1;
    while (v > 1 && static_cast<uint64_t>(v - 1) * static_cast<uint64_t>(v - 2) >= target) v--;
    while (static_cast<uint64_t>(v) * static_cast<uint64_t>(v - 1) < target) v++;
    return static_cast<uint64_t>(v);
}

static uint64_t plink_parallel_row_start(uint64_t n_ids, int parallel_idx /* 0-based */,
                                         int n_shards) {
    uint64_t ct_tot = n_ids * (n_ids - 1);  // "doubled" off-diagonal pair count
    uint64_t target = (ct_tot * static_cast<uint64_t>(parallel_idx)) / static_cast<uint64_t>(n_shards);
    uint64_t v = triangle_divide_off_diag(target);
    return v == 1 ? 0 : v;
}

// Seed for the shard-file-size verification below: plink_parallel_row_start
// is the real formula, but we still verify/search nearby in case of a plink
// version difference (or GCTA instead of plink) rather than trust it blindly.
static uint64_t estimate_row_start(uint64_t n_ids, int k, int n_shards) {
    return plink_parallel_row_start(n_ids, k - 1, n_shards);
}

static bool try_row_start(uint64_t candidate_start, uint64_t n_ids, uint64_t shard_floats,
                          uint64_t& row_end_out) {
    uint64_t consumed = 0;
    uint64_t i = candidate_start;
    while (consumed < shard_floats && i < n_ids) {
        consumed += (i + 1);
        i++;
    }
    if (consumed == shard_floats) {
        row_end_out = i;
        return true;
    }
    return false;
}

static ShardRange resolve_shard_range(uint64_t n_ids, int k, int n_shards, uint64_t shard_floats) {
    uint64_t guess = estimate_row_start(n_ids, k, n_shards);

    // Search a small window around the guess: covers both our formula being
    // an approximation and plink's own balancing not matching it exactly.
    static const int offsets[] = {0, -1, 1, -2, 2, -3, 3};
    for (int off : offsets) {
        int64_t candidate = static_cast<int64_t>(guess) + off;
        if (candidate < 0 || static_cast<uint64_t>(candidate) >= n_ids) continue;
        uint64_t row_end;
        if (try_row_start(static_cast<uint64_t>(candidate), n_ids, shard_floats, row_end)) {
            return {static_cast<uint64_t>(candidate), row_end};
        }
    }

    std::ostringstream oss;
    oss << "Could not determine this shard's row range from its file size.\n"
        << "  n_ids=" << n_ids << " parallel=" << k << "/" << n_shards
        << " shard_floats=" << shard_floats << " (best guess row_start=" << guess << ")\n"
        << "This usually means the shard file doesn't hold a whole number of "
           "complete rows, or the (k, n_shards) balancing plink used differs "
           "from ours by more than a few rows.";
    throw std::runtime_error(oss.str());
}

// ---------------------------------------------------------------------------
// accumulate
// ---------------------------------------------------------------------------

struct AccumulateArgs {
    std::string grm_id, shard, pheno, bins, out;
    int parallel_k = 0, parallel_n = 0, nblocks = 0;
    unsigned int seed = 1;
};

static AccumulateArgs parse_accumulate_args(int argc, char** argv) {
    AccumulateArgs a;
    for (int i = 2; i < argc; ++i) {
        std::string key = argv[i];
        auto need = [&](const std::string& flag) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("Missing value for " + flag);
            return argv[++i];
        };
        if (key == "--grm-id") a.grm_id = need(key);
        else if (key == "--shard") a.shard = need(key);
        else if (key == "--pheno") a.pheno = need(key);
        else if (key == "--bins") a.bins = need(key);
        else if (key == "--out") a.out = need(key);
        else if (key == "--nblocks") a.nblocks = std::stoi(need(key));
        else if (key == "--seed") a.seed = static_cast<unsigned int>(std::stoul(need(key)));
        else if (key == "--parallel") {
            a.parallel_k = std::stoi(need(key));
            a.parallel_n = std::stoi(need(key));
        } else throw std::runtime_error("Unknown argument: " + key);
    }
    if (a.grm_id.empty() || a.shard.empty() || a.pheno.empty() || a.bins.empty() ||
        a.out.empty() || a.parallel_k == 0 || a.parallel_n == 0 || a.nblocks <= 0) {
        throw std::runtime_error(
            "Usage: grm_shard_tool accumulate --grm-id FILE --shard FILE "
            "--parallel K N --pheno FILE --bins FILE --nblocks INT [--seed INT] --out FILE");
    }
    return a;
}

static int run_accumulate(int argc, char** argv) {
    AccumulateArgs args = parse_accumulate_args(argc, argv);

    std::vector<ID> ids = read_ids(args.grm_id);
    const uint64_t n_ids = ids.size();
    std::vector<double> y = read_pheno_aligned(args.pheno, ids);
    std::vector<Bin> bins = read_bins(args.bins);
    const int nbins = static_cast<int>(bins.size());
    std::vector<int> block_of = make_random_blocks(n_ids, args.nblocks, args.seed);

    std::ifstream shard(args.shard, std::ios::binary | std::ios::ate);
    if (!shard) throw std::runtime_error("Could not open shard file: " + args.shard);
    std::streamsize byte_size = shard.tellg();
    if (byte_size % 4 != 0) {
        throw std::runtime_error("Shard file size is not a multiple of 4 bytes: " + args.shard);
    }
    uint64_t shard_floats = static_cast<uint64_t>(byte_size) / 4;
    shard.seekg(0);

    ShardRange range = resolve_shard_range(n_ids, args.parallel_k, args.parallel_n, shard_floats);
    std::cerr << "[INFO] Shard " << args.parallel_k << "/" << args.parallel_n
              << " covers rows [" << range.row_start << ", " << range.row_end << ")\n";

    std::vector<Acc> full(nbins);
    std::vector<std::vector<Acc>> drop(args.nblocks, std::vector<Acc>(nbins));

    float g_f = 0.0f;
    for (uint64_t i = range.row_start; i < range.row_end; ++i) {
        for (uint64_t j = 0; j <= i; ++j) {
            shard.read(reinterpret_cast<char*>(&g_f), sizeof(float));
            if (!shard) throw std::runtime_error("Unexpected end of shard while reading row " +
                                                  std::to_string(i));
            if (i == j) continue;

            const double yi = y[i], yj = y[j];
            if (std::isnan(yi) || std::isnan(yj)) continue;

            const int k = get_bin(static_cast<double>(g_f), bins);
            if (k < 0) continue;

            const double prod = yi * yj;
            full[k].add(prod);

            const int bi = block_of[i], bj = block_of[j];
            drop[bi][k].add(prod);
            if (bj != bi) drop[bj][k].add(prod);
        }
    }

    std::ofstream out(args.out);
    if (!out) throw std::runtime_error("Could not open output file: " + args.out);
    out << "scope\tblock\tbin_index\tsum\tsum_sq\tn\n";
    out.precision(17);
    for (int k = 0; k < nbins; ++k) {
        out << "full\t-1\t" << k << '\t' << full[k].sum << '\t' << full[k].sum_sq << '\t'
            << full[k].n << '\n';
    }
    for (int b = 0; b < args.nblocks; ++b) {
        for (int k = 0; k < nbins; ++k) {
            if (drop[b][k].n == 0) continue;
            out << "drop\t" << b << '\t' << k << '\t' << drop[b][k].sum << '\t'
                << drop[b][k].sum_sq << '\t' << drop[b][k].n << '\n';
        }
    }
    std::cerr << "[INFO] Wrote accumulator to " << args.out << "\n";
    return 0;
}

// ---------------------------------------------------------------------------
// merge
// ---------------------------------------------------------------------------

struct MergeArgs {
    std::string acc_list, bins, out_prefix;
    int nblocks = 0;
};

static MergeArgs parse_merge_args(int argc, char** argv) {
    MergeArgs a;
    for (int i = 2; i < argc; ++i) {
        std::string key = argv[i];
        auto need = [&](const std::string& flag) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("Missing value for " + flag);
            return argv[++i];
        };
        if (key == "--acc-list") a.acc_list = need(key);
        else if (key == "--bins") a.bins = need(key);
        else if (key == "--out-prefix") a.out_prefix = need(key);
        else if (key == "--nblocks") a.nblocks = std::stoi(need(key));
        else throw std::runtime_error("Unknown argument: " + key);
    }
    if (a.acc_list.empty() || a.bins.empty() || a.out_prefix.empty() || a.nblocks <= 0) {
        throw std::runtime_error(
            "Usage: grm_shard_tool merge --acc-list FILE --bins FILE --nblocks INT --out-prefix PREFIX");
    }
    return a;
}

static void jackknife_summary_for_bin(int k, const std::vector<Acc>& full,
                                      const std::vector<std::vector<Acc>>& drop, int nblocks,
                                      double& jk_mean, double& jk_var, double& jk_se) {
    std::vector<double> means;
    means.reserve(nblocks);
    for (int b = 0; b < nblocks; ++b) {
        const double jk_sum = full[k].sum - drop[b][k].sum;
        const uint64_t jk_n = full[k].n - drop[b][k].n;
        if (jk_n == 0) {
            jk_mean = jk_var = jk_se = std::numeric_limits<double>::quiet_NaN();
            return;
        }
        means.push_back(jk_sum / static_cast<double>(jk_n));
    }

    double mean_theta = 0.0;
    for (double x : means) mean_theta += x;
    mean_theta /= static_cast<double>(nblocks);

    double ssd = 0.0;
    for (double x : means) {
        const double d = x - mean_theta;
        ssd += d * d;
    }

    jk_mean = mean_theta;
    jk_var = ((static_cast<double>(nblocks) - 1.0) / static_cast<double>(nblocks)) * ssd;
    jk_se = std::sqrt(jk_var);
}

static int run_merge(int argc, char** argv) {
    MergeArgs args = parse_merge_args(argc, argv);
    std::vector<Bin> bins = read_bins(args.bins);
    const int nbins = static_cast<int>(bins.size());

    std::vector<Acc> full(nbins);
    std::vector<std::vector<Acc>> drop(args.nblocks, std::vector<Acc>(nbins));

    std::ifstream list(args.acc_list);
    if (!list) throw std::runtime_error("Could not open shard list file: " + args.acc_list);

    std::string acc_path;
    int n_shards_read = 0;
    while (std::getline(list, acc_path)) {
        if (acc_path.empty()) continue;
        std::ifstream in(acc_path);
        if (!in) throw std::runtime_error("Could not open accumulator file: " + acc_path);

        std::string line;
        std::getline(in, line);  // header
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            std::istringstream iss(line);
            std::string scope;
            int block, bin_index;
            double sum, sum_sq;
            uint64_t n;
            iss >> scope >> block >> bin_index >> sum >> sum_sq >> n;

            if (bin_index < 0 || bin_index >= nbins) {
                throw std::runtime_error("Bin index out of range in " + acc_path +
                                         " -- bins file mismatch between accumulate and merge?");
            }
            Acc piece{sum, sum_sq, n};
            if (scope == "full") {
                full[bin_index].add(piece);
            } else if (scope == "drop") {
                if (block < 0 || block >= args.nblocks) {
                    throw std::runtime_error("Block index out of range in " + acc_path +
                                             " -- nblocks mismatch between accumulate and merge?");
                }
                drop[block][bin_index].add(piece);
            }
        }
        n_shards_read++;
    }
    std::cerr << "[INFO] Merged " << n_shards_read << " shard accumulator files\n";

    std::ofstream out_full(args.out_prefix + ".full.tsv");
    std::ofstream out_jk(args.out_prefix + ".jk.tsv");
    out_full << "bin_index\tbin_left\tbin_right\tbin_midpoint\tfull_sum\tfull_sum_sq\tfull_n\t"
                "full_mean\tfull_sd\tfull_se\tjk_mean\tjk_var\tjk_se\n";
    out_jk << "block\tbin_index\tbin_left\tbin_right\tbin_midpoint\tjk_sum\tjk_sum_sq\tjk_n\t"
              "jk_mean\tjk_sd\tjk_se\n";

    for (int k = 0; k < nbins; ++k) {
        const double mean = safe_mean(full[k].sum, full[k].n);
        const double sd = safe_sample_sd(full[k].sum, full[k].sum_sq, full[k].n);
        const double se = safe_se_of_mean(full[k].sum, full[k].sum_sq, full[k].n);

        double jk_mean = std::numeric_limits<double>::quiet_NaN();
        double jk_var = std::numeric_limits<double>::quiet_NaN();
        double jk_se = std::numeric_limits<double>::quiet_NaN();
        jackknife_summary_for_bin(k, full, drop, args.nblocks, jk_mean, jk_var, jk_se);

        out_full << (k + 1) << '\t' << bins[k].left << '\t' << bins[k].right << '\t'
                 << bin_midpoint(bins[k]) << '\t' << full[k].sum << '\t' << full[k].sum_sq << '\t'
                 << full[k].n << '\t' << mean << '\t' << sd << '\t' << se << '\t' << jk_mean << '\t'
                 << jk_var << '\t' << jk_se << '\n';
    }

    for (int b = 0; b < args.nblocks; ++b) {
        for (int k = 0; k < nbins; ++k) {
            const double jk_sum = full[k].sum - drop[b][k].sum;
            const double jk_sum_sq = full[k].sum_sq - drop[b][k].sum_sq;
            const uint64_t jk_n = full[k].n - drop[b][k].n;

            const double jk_mean = safe_mean(jk_sum, jk_n);
            const double jk_sd = safe_sample_sd(jk_sum, jk_sum_sq, jk_n);
            const double jk_se = safe_se_of_mean(jk_sum, jk_sum_sq, jk_n);

            out_jk << b << '\t' << (k + 1) << '\t' << bins[k].left << '\t' << bins[k].right << '\t'
                   << bin_midpoint(bins[k]) << '\t' << jk_sum << '\t' << jk_sum_sq << '\t' << jk_n
                   << '\t' << jk_mean << '\t' << jk_sd << '\t' << jk_se << '\n';
        }
    }

    std::cerr << "[INFO] Wrote " << args.out_prefix << ".full.tsv and " << args.out_prefix
              << ".jk.tsv\n";
    return 0;
}

// ---------------------------------------------------------------------------
// ranges
// ---------------------------------------------------------------------------

static int run_ranges(int argc, char** argv) {
    uint64_t n_ids = 0;
    int k = 0, n_shards = 0;
    for (int i = 2; i < argc; ++i) {
        std::string key = argv[i];
        auto need = [&](const std::string& flag) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("Missing value for " + flag);
            return argv[++i];
        };
        if (key == "--n-ids") n_ids = std::stoull(need(key));
        else if (key == "--parallel") { k = std::stoi(need(key)); n_shards = std::stoi(need(key)); }
        else throw std::runtime_error("Unknown argument: " + key);
    }
    if (n_ids == 0 || k == 0 || n_shards == 0) {
        throw std::runtime_error("Usage: grm_shard_tool ranges --n-ids INT --parallel K N");
    }

    uint64_t row_start = estimate_row_start(n_ids, k, n_shards);
    uint64_t row_end = estimate_row_start(n_ids, k + 1, n_shards);
    if (k == n_shards) row_end = n_ids;
    uint64_t n_floats = cumulative_entries(row_end) - cumulative_entries(row_start);

    std::cout << "shard " << k << "/" << n_shards << ": rows [" << row_start << ", " << row_end
              << "), " << n_floats << " floats (" << n_floats * 4 << " bytes), "
              << "byte offset of first row = " << cumulative_entries(row_start) * 4 << "\n";
    return 0;
}

// ---------------------------------------------------------------------------

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: grm_shard_tool <accumulate|merge|ranges> [options]\n";
        return 1;
    }
    std::string cmd = argv[1];
    try {
        if (cmd == "accumulate") return run_accumulate(argc, argv);
        if (cmd == "merge") return run_merge(argc, argv);
        if (cmd == "ranges") return run_ranges(argc, argv);
        std::cerr << "Unknown subcommand: " << cmd << "\n";
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}
