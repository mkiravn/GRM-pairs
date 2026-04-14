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
#include <unordered_set>
#include <utility>
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
};

struct Args {
    std::string grm_prefix;
    std::string pheno_path;
    std::string bins_path;
    std::string out_prefix;
    std::string exclude_pairs_path;
    int nblocks = 100;
    unsigned int seed = 1;
};

static Args parse_args(int argc, char** argv) {
    Args args;

    for (int i = 1; i < argc; ++i) {
        std::string key = argv[i];

        auto need_value = [&](const std::string& flag) -> std::string {
            if (i + 1 >= argc) {
                throw std::runtime_error("Missing value for " + flag);
            }
            return argv[++i];
        };

        if (key == "--grm-prefix") {
            args.grm_prefix = need_value(key);
        } else if (key == "--pheno") {
            args.pheno_path = need_value(key);
        } else if (key == "--bins") {
            args.bins_path = need_value(key);
        } else if (key == "--out-prefix") {
            args.out_prefix = need_value(key);
        } else if (key == "--exclude-pairs") {
            args.exclude_pairs_path = need_value(key);
        } else if (key == "--nblocks") {
            args.nblocks = std::stoi(need_value(key));
        } else if (key == "--seed") {
            args.seed = static_cast<unsigned int>(std::stoul(need_value(key)));
        } else {
            throw std::runtime_error("Unknown argument: " + key);
        }
    }

    if (args.grm_prefix.empty() ||
        args.pheno_path.empty() ||
        args.bins_path.empty() ||
        args.out_prefix.empty()) {
        throw std::runtime_error(
            "Usage: --grm-prefix PREFIX --pheno FILE --bins FILE "
            "--out-prefix PREFIX [--exclude-pairs FILE] [--nblocks 100] [--seed 1]"
        );
    }

    if (args.nblocks <= 0) {
        throw std::runtime_error("--nblocks must be > 0");
    }

    return args;
}

static bool is_missing_token(const std::string& s) {
    return s == "NA" || s == "NaN" || s == "nan" || s == ".";
}

static double parse_double_or_nan(const std::string& s) {
    if (is_missing_token(s)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    try {
        size_t idx = 0;
        double x = std::stod(s, &idx);
        if (idx != s.size()) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return x;
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

static std::vector<ID> read_ids(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Could not open grm.id file: " + path);
    }

    std::vector<ID> ids;
    std::string fid, iid;
    while (in >> fid >> iid) {
        ids.push_back({fid, iid});
    }

    if (ids.empty()) {
        throw std::runtime_error("No IDs found in grm.id file: " + path);
    }

    return ids;
}

static std::vector<double> read_pheno_aligned(
    const std::string& path,
    const std::vector<ID>& ids
) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Could not open phenotype file: " + path);
    }

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
        if (it != pheno.end()) {
            y[i] = it->second;
        }
    }

    return y;
}

static std::vector<Bin> read_bins(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Could not open bins file: " + path);
    }

    std::vector<Bin> bins;
    double a, b;
    while (in >> a >> b) {
        if (!(a < b)) {
            throw std::runtime_error("Invalid bin with left >= right");
        }
        bins.push_back({a, b});
    }

    if (bins.empty()) {
        throw std::runtime_error("No bins found in file: " + path);
    }

    std::sort(bins.begin(), bins.end(), [](const Bin& x, const Bin& y) {
        return x.left < y.left;
    });

    return bins;
}

static int get_bin(double g, const std::vector<Bin>& bins) {
    int lo = 0;
    int hi = static_cast<int>(bins.size()) - 1;

    while (lo <= hi) {
        int mid = lo + (hi - lo) / 2;
        const Bin& b = bins[mid];

        bool in_bin = false;
        if (mid == static_cast<int>(bins.size()) - 1) {
            in_bin = (g >= b.left && g <= b.right);
        } else {
            in_bin = (g >= b.left && g < b.right);
        }

        if (in_bin) return mid;
        if (g < b.left) {
            hi = mid - 1;
        } else {
            lo = mid + 1;
        }
    }

    return -1;
}

static std::vector<int> make_random_blocks(std::size_t n, int nblocks, unsigned int seed) {
    if (static_cast<std::size_t>(nblocks) > n) {
        throw std::runtime_error("Number of blocks exceeds number of individuals");
    }

    std::vector<std::size_t> perm(n);
    for (std::size_t i = 0; i < n; ++i) {
        perm[i] = i;
    }

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

static double bin_midpoint(const Bin& b) {
    return 0.5 * (b.left + b.right);
}

static double safe_mean(double sum, std::uint64_t n) {
    if (n == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return sum / static_cast<double>(n);
}

static double safe_sample_variance(double sum, double sum_sq, std::uint64_t n) {
    if (n < 2) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double nd = static_cast<double>(n);
    const double numer = sum_sq - (sum * sum) / nd;
    double var = numer / static_cast<double>(n - 1);

    if (var < 0.0 && std::abs(var) < 1e-12) {
        var = 0.0;
    }
    return var;
}

static double safe_sample_sd(double sum, double sum_sq, std::uint64_t n) {
    const double var = safe_sample_variance(sum, sum_sq, n);
    if (std::isnan(var)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::sqrt(var);
}

static double safe_se_of_mean(double sum, double sum_sq, std::uint64_t n) {
    const double var = safe_sample_variance(sum, sum_sq, n);
    if (std::isnan(var)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::sqrt(var / static_cast<double>(n));
}

static void jackknife_summary_for_bin(
    std::size_t k,
    const std::vector<Acc>& full,
    const std::vector<std::vector<Acc>>& excluded,
    int nblocks,
    double& jk_mean,
    double& jk_var,
    double& jk_se
) {
    std::vector<double> means;
    means.reserve(nblocks);

    for (int b = 0; b < nblocks; ++b) {
        const double jk_sum = full[k].sum - excluded[b][k].sum;
        const std::uint64_t jk_n = full[k].n - excluded[b][k].n;
        if (jk_n == 0) {
            jk_mean = std::numeric_limits<double>::quiet_NaN();
            jk_var = std::numeric_limits<double>::quiet_NaN();
            jk_se = std::numeric_limits<double>::quiet_NaN();
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

static void write_full_output(
    const std::string& path,
    const std::vector<Bin>& bins,
    const std::vector<Acc>& full,
    const std::vector<std::vector<Acc>>& excluded,
    int nblocks
) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Could not open full output file: " + path);
    }

    out << "bin_index\tbin_left\tbin_right\tbin_midpoint\tfull_sum\tfull_sum_sq\tfull_n\tfull_mean\tfull_sd\tfull_se\tjk_mean\tjk_var\tjk_se\n";

    for (std::size_t k = 0; k < bins.size(); ++k) {
        const double mean = safe_mean(full[k].sum, full[k].n);
        const double sd   = safe_sample_sd(full[k].sum, full[k].sum_sq, full[k].n);
        const double se   = safe_se_of_mean(full[k].sum, full[k].sum_sq, full[k].n);

        double jk_mean = std::numeric_limits<double>::quiet_NaN();
        double jk_var  = std::numeric_limits<double>::quiet_NaN();
        double jk_se   = std::numeric_limits<double>::quiet_NaN();

        jackknife_summary_for_bin(k, full, excluded, nblocks, jk_mean, jk_var, jk_se);

        out << (k + 1) << '\t'
            << bins[k].left << '\t'
            << bins[k].right << '\t'
            << bin_midpoint(bins[k]) << '\t'
            << full[k].sum << '\t'
            << full[k].sum_sq << '\t'
            << full[k].n << '\t'
            << mean << '\t'
            << sd << '\t'
            << se << '\t'
            << jk_mean << '\t'
            << jk_var << '\t'
            << jk_se << '\n';
    }
}

static void write_jk_output(
    const std::string& path,
    const std::vector<Bin>& bins,
    const std::vector<Acc>& full,
    const std::vector<std::vector<Acc>>& excluded,
    int nblocks
) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Could not open jackknife output file: " + path);
    }

    out << "block\tbin_index\tbin_left\tbin_right\tbin_midpoint\tjk_sum\tjk_sum_sq\tjk_n\tjk_mean\tjk_sd\tjk_se\n";

    for (int b = 0; b < nblocks; ++b) {
        for (std::size_t k = 0; k < bins.size(); ++k) {
            const double jk_sum = full[k].sum - excluded[b][k].sum;
            const double jk_sum_sq = full[k].sum_sq - excluded[b][k].sum_sq;
            const std::uint64_t jk_n = full[k].n - excluded[b][k].n;

            const double jk_mean = safe_mean(jk_sum, jk_n);
            const double jk_sd   = safe_sample_sd(jk_sum, jk_sum_sq, jk_n);
            const double jk_se   = safe_se_of_mean(jk_sum, jk_sum_sq, jk_n);

            out << b << '\t'
                << (k + 1) << '\t'
                << bins[k].left << '\t'
                << bins[k].right << '\t'
                << bin_midpoint(bins[k]) << '\t'
                << jk_sum << '\t'
                << jk_sum_sq << '\t'
                << jk_n << '\t'
                << jk_mean << '\t'
                << jk_sd << '\t'
                << jk_se << '\n';
        }
    }
}

struct PairKey {
    std::string a;
    std::string b;

    bool operator==(const PairKey& other) const {
        return a == other.a && b == other.b;
    }
};

struct PairKeyHash {
    std::size_t operator()(const PairKey& x) const {
        return std::hash<std::string>()(x.a) ^ (std::hash<std::string>()(x.b) << 1);
    }
};

static PairKey make_pair_key(const std::string& x, const std::string& y) {
    if (x <= y) return {x, y};
    return {y, x};
}

static std::unordered_set<PairKey, PairKeyHash> read_excluded_pairs(const std::string& path) {
    std::unordered_set<PairKey, PairKeyHash> excluded;

    if (path.empty()) {
        return excluded;
    }

    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Could not open exclude-pairs file: " + path);
    }

    std::string line;
    bool first_line = true;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string x, y;

        if (!(iss >> x >> y)) continue;

        // skip likely header
        if (first_line) {
            first_line = false;
            if ((x == "eid1" && y == "eid2") ||
                (x == "IID1" && y == "IID2") ||
                (x == "participant.eid1" && y == "participant.eid2")) {
                continue;
            }
        }

        excluded.insert(make_pair_key(x, y));
    }

    return excluded;
}

int main(int argc, char** argv) {
    try {
        Args args = parse_args(argc, argv);

        const std::string id_path = args.grm_prefix + ".grm.id";
        const std::string grm_path = args.grm_prefix + ".grm.bin";
        const std::string out_full = args.out_prefix + ".full.tsv";
        const std::string out_jk = args.out_prefix + ".jk.tsv";

        std::cerr << "[INFO] Reading IDs from " << id_path << "\n";
        std::vector<ID> ids = read_ids(id_path);
        const std::size_t n = ids.size();

        std::cerr << "[INFO] Number of individuals: " << n << "\n";

        std::cerr << "[INFO] Reading phenotype and aligning to grm.id order\n";
        std::vector<double> y = read_pheno_aligned(args.pheno_path, ids);

        std::cerr << "[INFO] Reading bins from " << args.bins_path << "\n";
        std::vector<Bin> bins = read_bins(args.bins_path);
        const int nbins = static_cast<int>(bins.size());

        std::cerr << "[INFO] Assigning individuals randomly to " << args.nblocks
                  << " blocks with seed " << args.seed << "\n";
        std::vector<int> block_of = make_random_blocks(n, args.nblocks, args.seed);

        std::cerr << "[INFO] Reading excluded pairs\n";
        std::unordered_set<PairKey, PairKeyHash> excluded_pairs =
            read_excluded_pairs(args.exclude_pairs_path);
        std::cerr << "[INFO] Excluded pair count: " << excluded_pairs.size() << "\n";

        std::vector<Acc> full(nbins);
        std::vector<std::vector<Acc>> excluded(
            args.nblocks, std::vector<Acc>(nbins)
        );

        std::cerr << "[INFO] Streaming GRM from " << grm_path << "\n";
        std::ifstream grm(grm_path, std::ios::binary);
        if (!grm) {
            throw std::runtime_error("Could not open grm.bin file: " + grm_path);
        }

        float g_f = 0.0f;

        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j <= i; ++j) {
                grm.read(reinterpret_cast<char*>(&g_f), sizeof(float));
                if (!grm) {
                    throw std::runtime_error("Unexpected end of grm.bin");
                }

                if (i == j) {
                    continue;
                }

                // Skip explicitly excluded pairs, matching on IID only
                if (!excluded_pairs.empty()) {
                    PairKey key = make_pair_key(ids[i].iid, ids[j].iid);
                    if (excluded_pairs.find(key) != excluded_pairs.end()) {
                        continue;
                    }
                }

                const double g = static_cast<double>(g_f);
                const int k = get_bin(g, bins);
                if (k < 0) {
                    continue;
                }

                const double yi = y[i];
                const double yj = y[j];
                if (std::isnan(yi) || std::isnan(yj)) {
                    continue;
                }

                const double prod = yi * yj;
                const double prod_sq = prod * prod;

                full[k].sum += prod;
                full[k].sum_sq += prod_sq;
                full[k].n += 1;

                const int bi = block_of[i];
                const int bj = block_of[j];

                excluded[bi][k].sum += prod;
                excluded[bi][k].sum_sq += prod_sq;
                excluded[bi][k].n += 1;

                if (bj != bi) {
                    excluded[bj][k].sum += prod;
                    excluded[bj][k].sum_sq += prod_sq;
                    excluded[bj][k].n += 1;
                }
            }
        }

        std::cerr << "[INFO] Writing full output to " << out_full << "\n";
        write_full_output(out_full, bins, full, excluded, args.nblocks);

        std::cerr << "[INFO] Writing jackknife output to " << out_jk << "\n";
        write_jk_output(out_jk, bins, full, excluded, args.nblocks);

        std::cerr << "[INFO] Done\n";
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "[ERROR] " << e.what() << "\n";
        return 1;
    }
}