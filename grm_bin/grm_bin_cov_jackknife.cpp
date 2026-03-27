#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct BinTotals {
    uint64_t count = 0;
    double sum_grm = 0.0;
    double sum_cov = 0.0;
};

static std::vector<double> make_bins(double min_edge = -0.3,
                                     double split1 = -0.02,
                                     double split2 = 0.02,
                                     double max_edge = 1.05,
                                     double small_width = 0.001,
                                     double large_width = 0.005) {
    std::vector<double> edges;
    edges.push_back(min_edge);
    edges.push_back(split1);

    for (double x = split1 + small_width; x <= split2 + 1e-12; x += small_width) {
        edges.push_back(x);
    }
    for (double x = split2 + large_width; x <= max_edge + 1e-12; x += large_width) {
        edges.push_back(x);
    }
    if (std::abs(edges.back() - max_edge) > 1e-12) {
        edges.push_back(max_edge);
    }
    return edges;
}

static int find_bin(double x, const std::vector<double>& edges) {
    if (x < edges.front() || x > edges.back()) return -1;
    if (x == edges.back()) return static_cast<int>(edges.size()) - 2;
    auto it = std::upper_bound(edges.begin(), edges.end(), x);
    int idx = static_cast<int>(it - edges.begin()) - 1;
    if (idx < 0 || idx >= static_cast<int>(edges.size()) - 1) return -1;
    return idx;
}

static size_t count_grm_ids(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) {
        throw std::runtime_error("Could not open GRM ID file: " + filename);
    }

    size_t n = 0;
    std::string fid, iid;
    while (in >> fid >> iid) {
        ++n;
    }

    if (n == 0) {
        throw std::runtime_error("No IDs found in GRM ID file: " + filename);
    }

    return n;
}

static bool parse_line(const std::string& line, double& y) {
    if (line.empty()) return false;

    std::istringstream iss(line);
    std::string id, token;
    if (!(iss >> id >> token)) {
        return false;
    }

    if (token == "NA" || token == "NaN" || token == "nan" || token == "-999") {
        return false;
    }

    try {
        y = std::stod(token);
        return true;
    } catch (...) {
        return false;
    }
}

static std::vector<int> assign_blocks(size_t n, int n_blocks) {
    if (n_blocks < 2) {
        throw std::runtime_error("n_blocks must be at least 2.");
    }

    std::vector<int> block(n, 0);
    for (size_t i = 0; i < n; ++i) {
        block[i] = static_cast<int>((static_cast<long double>(i) * n_blocks) / n);
        if (block[i] >= n_blocks) block[i] = n_blocks - 1;
    }
    return block;
}

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <aligned_pheno.txt> <grm_prefix> <output.tsv> <n_blocks> [jackknife_detail.tsv]\n";
        std::cerr << "aligned_pheno.txt should have two columns: ID phenotype, in GRM order.\n";
        std::cerr << "Missing phenotypes may be encoded as NA or -999.\n";
        return 1;
    }

    const std::string pheno_file = argv[1];
    const std::string grm_prefix = argv[2];
    const std::string output_file = argv[3];
    const int n_blocks = std::stoi(argv[4]);
    const bool write_detail = (argc >= 6);
    const std::string detail_file = write_detail ? argv[5] : "";

    const std::string grm_id_file = grm_prefix + ".grm.id";
    const std::string grm_bin_file = grm_prefix + ".grm.bin";

    try {
        const size_t n = count_grm_ids(grm_id_file);

        std::cerr << "Phenotype file: " << pheno_file << "\n";
        std::cerr << "GRM prefix:     " << grm_prefix << "\n";
        std::cerr << "n individuals:  " << n << "\n";
        std::cerr << "n blocks:       " << n_blocks << "\n";

        std::vector<double> y(n, 0.0);
        std::vector<uint8_t> keep(n, 0);

        {
            std::ifstream in(pheno_file);
            if (!in) {
                throw std::runtime_error("Could not open phenotype file: " + pheno_file);
            }

            std::string line;
            size_t i = 0;
            size_t nonmissing = 0;

            while (std::getline(in, line)) {
                if (line.empty()) continue;

                if (i >= n) {
                    throw std::runtime_error("Phenotype file has more rows than .grm.id");
                }

                double val = 0.0;
                if (parse_line(line, val)) {
                    y[i] = val;
                    keep[i] = 1;
                    ++nonmissing;
                }

                ++i;
            }

            if (i != n) {
                std::ostringstream oss;
                oss << "Phenotype file has " << i
                    << " rows, but .grm.id has " << n << " rows.";
                throw std::runtime_error(oss.str());
            }

            if (nonmissing < 2) {
                throw std::runtime_error("Need at least 2 non-missing phenotypes.");
            }

            std::cerr << "Non-missing phenotypes: " << nonmissing << "\n";
            std::cerr << "Phenotype assumed pre-standardized/residualized; values left unchanged.\n";
        }

        const auto edges = make_bins();
        const size_t nbins = edges.size() - 1;

        std::vector<BinTotals> full(nbins);
        std::vector<std::vector<BinTotals>> drop(
            n_blocks, std::vector<BinTotals>(nbins)
        );

        const std::vector<int> block = assign_blocks(n, n_blocks);

        std::cerr << "Number of bins: " << nbins << "\n";
        std::cerr << "Range: [" << edges.front() << ", " << edges.back() << "]\n";

        std::ifstream grm(grm_bin_file, std::ios::binary);
        if (!grm) {
            throw std::runtime_error("Could not open GRM binary file: " + grm_bin_file);
        }

        uint64_t total_pairs_seen = 0;
        uint64_t total_pairs_used = 0;
        uint64_t out_of_range = 0;

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j <= i; ++j) {
                float grm_val_f;
                grm.read(reinterpret_cast<char*>(&grm_val_f), sizeof(float));

                if (!grm) {
                    throw std::runtime_error("Unexpected EOF while reading " + grm_bin_file);
                }

                ++total_pairs_seen;

                if (i == j) continue;
                if (!keep[i] || !keep[j]) continue;

                const double grm_val = static_cast<double>(grm_val_f);
                const int k = find_bin(grm_val, edges);

                if (k < 0) {
                    ++out_of_range;
                    continue;
                }

                const double cp = y[i] * y[j];
                const int bi = block[i];
                const int bj = block[j];

                full[k].count += 1;
                full[k].sum_grm += grm_val;
                full[k].sum_cov += cp;

                drop[bi][k].count += 1;
                drop[bi][k].sum_grm += grm_val;
                drop[bi][k].sum_cov += cp;

                if (bj != bi) {
                    drop[bj][k].count += 1;
                    drop[bj][k].sum_grm += grm_val;
                    drop[bj][k].sum_cov += cp;
                }

                ++total_pairs_used;
            }
        }

        char extra;
        if (grm.read(&extra, 1)) {
            std::cerr << "Warning: extra bytes found after expected end of GRM file.\n";
        }

        std::cerr << "Total GRM entries seen (including diagonal): "
                  << total_pairs_seen << "\n";
        std::cerr << "Pairs used (off-diagonal, non-missing phenotype): "
                  << total_pairs_used << "\n";

        if (out_of_range > 0) {
            std::cerr << "Warning: " << out_of_range
                      << " pairs had GRM values outside bin range and were skipped.\n";
        }

        std::vector<double> full_avg_grm(nbins, std::numeric_limits<double>::quiet_NaN());
        std::vector<double> full_avg_cov(nbins, std::numeric_limits<double>::quiet_NaN());
        std::vector<double> se_cov(nbins, std::numeric_limits<double>::quiet_NaN());
        std::vector<std::vector<double>> loo_avg_grm(
            n_blocks, std::vector<double>(nbins, std::numeric_limits<double>::quiet_NaN())
        );
        std::vector<std::vector<double>> loo_avg_cov(
            n_blocks, std::vector<double>(nbins, std::numeric_limits<double>::quiet_NaN())
        );
        std::vector<std::vector<uint64_t>> loo_count(
            n_blocks, std::vector<uint64_t>(nbins, 0)
        );

        for (size_t k = 0; k < nbins; ++k) {
            if (full[k].count == 0) continue;

            full_avg_grm[k] = full[k].sum_grm / static_cast<double>(full[k].count);
            full_avg_cov[k] = full[k].sum_cov / static_cast<double>(full[k].count);

            std::vector<double> loo_cov;
            loo_cov.reserve(n_blocks);

            for (int b = 0; b < n_blocks; ++b) {
                const uint64_t n_keep = full[k].count - drop[b][k].count;
                if (n_keep == 0) continue;

                const double grm_keep = full[k].sum_grm - drop[b][k].sum_grm;
                const double s_keep = full[k].sum_cov - drop[b][k].sum_cov;
                loo_avg_grm[b][k] = grm_keep / static_cast<double>(n_keep);
                loo_avg_cov[b][k] = s_keep / static_cast<double>(n_keep);
                loo_count[b][k] = n_keep;
                loo_cov.push_back(s_keep / static_cast<double>(n_keep));
            }

            if (loo_cov.size() < 2) continue;

            double mean_loo = 0.0;
            for (double v : loo_cov) mean_loo += v;
            mean_loo /= static_cast<double>(loo_cov.size());

            double var = 0.0;
            for (double v : loo_cov) {
                const double d = v - mean_loo;
                var += d * d;
            }

            var *= (static_cast<double>(loo_cov.size() - 1) /
                    static_cast<double>(loo_cov.size()));

            se_cov[k] = std::sqrt(var);
        }

        std::ofstream out(output_file);
        if (!out) {
            throw std::runtime_error("Could not open output file: " + output_file);
        }

        out << "bin\tlower\tupper\tavg_grm\tavg_pheno_crossprod\tn_pairs\tse\n";
        out << std::fixed << std::setprecision(6);

        for (size_t k = 0; k < nbins; ++k) {
            if (full[k].count == 0) continue;

            out << (k + 1) << '\t'
                << edges[k] << '\t'
                << edges[k + 1] << '\t'
                << full_avg_grm[k] << '\t'
                << full_avg_cov[k] << '\t'
                << full[k].count << '\t'
                << se_cov[k] << '\n';

            std::cout << std::setw(5) << (k + 1) << "  "
                      << std::setw(9) << edges[k] << "  "
                      << std::setw(9) << edges[k + 1] << "  "
                      << std::setw(12) << full_avg_grm[k] << "  "
                      << std::setw(14) << full_avg_cov[k] << "  "
                      << std::setw(12) << full[k].count << "  "
                      << std::setw(12) << se_cov[k] << "\n";
        }

        std::cerr << "Wrote output to " << output_file << "\n";

        if (write_detail) {
            std::ofstream detail(detail_file);
            if (!detail) {
                throw std::runtime_error("Could not open detail output file: " + detail_file);
            }

            detail << "block\tbin\tlower\tupper\tavg_grm\tavg_pheno_crossprod\tn_pairs\n";
            detail << std::fixed << std::setprecision(6);

            for (int b = 0; b < n_blocks; ++b) {
                for (size_t k = 0; k < nbins; ++k) {
                    if (loo_count[b][k] == 0) continue;
                    detail << (b + 1) << '\t'
                           << (k + 1) << '\t'
                           << edges[k] << '\t'
                           << edges[k + 1] << '\t'
                           << loo_avg_grm[b][k] << '\t'
                           << loo_avg_cov[b][k] << '\t'
                           << loo_count[b][k] << '\n';
                }
            }

            std::cerr << "Wrote jackknife detail to " << detail_file << "\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
