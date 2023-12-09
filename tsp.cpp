#pragma GCC optimize("O3,unroll-loops")
#include <bits/stdc++.h>
using namespace std::chrono;
using namespace std::literals::chrono_literals;

using index_t = uint_fast16_t;
using weight_t = std::int64_t;

template <typename T>
bool in_list(T x) { return 0; }

template <typename T, typename... Args>
bool in_list(T x, T a, Args... args) {
    return x == a || in_list(x, args...);
}

template <class T>
class Matrix {
    std::vector<T> data;
    size_t _m;
    size_t _n;

public:
    Matrix() = default;
    Matrix(size_t m, size_t n) : _m(m), _n(n), data(m * n) {}

    const T& operator()(size_t i, size_t j) const {
        return data[i * _n + j];
    }
    T& operator()(size_t i, size_t j) {
        return data[i * _n + j];
    }

    size_t m() const { return _m; }
    size_t n() const { return _n; }
};

struct Point {
    double x;
    double y;

    double d(const Point& other) const {
        double dx = x - other.x;
        double dy = y - other.y;
        return sqrt(dx * dx + dy * dy);
    }
};

// 2-level or array representation?
// Depends on rate of convergence
struct Tour {
    std::vector<index_t> t;
    std::vector<size_t> pos;

    // To avoid bugs, we reset dlb aggressively
    // Moreover dlb2 does not imply dlb3
    struct DLB {
        bool prev2 = 0;
        bool next2 = 0;
        bool prev3 = 0;
        bool next3 = 0;

        void reverse() {
            std::swap(prev2, next2);
            std::swap(prev3, next3);
        }

        void reset() {
            prev2 = next2 = 0;
            prev3 = next3 = 0;
        }

        bool try_skip_2() {
            if (prev2 && next2) return 1;
            next2 = 1;  // not prev
            return 0;
        }

        bool try_skip_3() {
            if (prev3 && next3) return 1;
            next3 = 1;  // not prev
            return 0;
        }
    };
    std::vector<DLB> dlb;

    Tour() = default;
    Tour(std::vector<index_t>&& v) : t(v) {
        size_t n = v.size();
        pos.resize(n);
        dlb.resize(n);
        for (size_t i = 0; i < n; ++i) pos[v[i]] = i;
    }

    size_t n() const {
        return t.size();
    }

    size_t inc(size_t x) const {
        ++x;
        if (x >= n()) return x - n();
        return x;
    }

    size_t dec(size_t x) const {
        x += n() - 1;
        if (x >= n()) return x - n();
        return x;
    }

    size_t segment_length(index_t t1, index_t t2) const {
        size_t start = pos[t1];
        size_t end = pos[t2];
        return start <= end ? end - start : (end + n()) - start;
    }

    // Copied reverse method from: https://github.com/estan/tsp/blob/master/tsp.cpp
    void reverse(index_t t1, index_t t2) {
        /*
            skip[prev(t1)] = 0;
            for (size_t i = pos[t2]; i != pos[t1]; i = dec(i)) {
                skip[t[i]] = skip[t[dec(i)]];
            }
        */
        size_t m = (segment_length(t1, t2) + 1) / 2;
        size_t i = pos[t1];
        size_t j = pos[t2];
        while (m--) { 
            std::swap(t[i], t[j]);
            dlb[t[i]].reverse();
            dlb[t[j]].reverse();
            pos[t[i]] = i;
            pos[t[j]] = j;
            i = inc(i);
            j = dec(j);
        }
    }

    index_t next(index_t x) const {
        return t[inc(pos[x])];
    }

    index_t next_i(size_t i) const {
        return t[inc(i)];
    }

    index_t prev(index_t x) const {
        return t[dec(pos[x])];
    }

    index_t prev_i(size_t i) const {
        return t[dec(i)];
    }

    void reset_dlb(index_t x) {
        dlb[x].reset();
        dlb[prev(x)].reset();
        dlb[next(x)].reset();
    }

    // Interface
    const std::vector<index_t>& vec() const {
        return t;
    }

    void perturb(std::mt19937& rng) {
        if (n() < 8) return;

        std::vector<index_t> v;
        v.reserve(n());

        // Segments [a, b)
        std::uniform_int_distribution<size_t> d0(1, n() / 4);
        size_t i1 = d0(rng);
        size_t i3 = d0(rng) + i1;
        size_t i5 = d0(rng) + i3;
        size_t i7 = n() - 1;

        // We reset more than necessary here...
        for (auto i : {i1, i3, i5, i7}) reset_dlb(t[i]);

        std::copy(t.begin(),       t.begin() + i1,    back_inserter(v));
        std::copy(t.begin() + i5,  t.end(),           back_inserter(v));
        std::copy(t.begin() + i3,  t.begin() + i5,    back_inserter(v));
        std::copy(t.begin() + i1,  t.begin() + i3,    back_inserter(v));

        if (n() <= 30) {
            std::shuffle(v.begin(), v.end(), rng);
        } else {
            std::uniform_int_distribution<size_t> d1(0, n() - 30);
            size_t l = d1(rng);
            auto is = v.begin() + l;
            auto ie = is + 30;
            std::shuffle(is, ie, rng);
            for (auto it = is; it != ie; ++it) reset_dlb(*it);
        }

        t = std::move(v);
        for (size_t i = 0; i < n(); ++i) pos[t[i]] = i;
    }

    // 12a - 34 -> 1a - 324
    void swap2b(
        index_t t1, index_t t2,
        index_t t3, index_t t4,
        index_t b
    ) {
        for (size_t i = pos[t2]; inc(i) != pos[t4]; i = inc(i)) {
            std::swap(pos[t[i]], pos[t[inc(i)]]);
            std::swap(t[i], t[inc(i)]);
        }
    }

    // 12 - a34 -> 132 - a4
    void swap2c(
        index_t t1, index_t t2,
        index_t t3, index_t t4,
        index_t c
    ) {
        for (size_t i = pos[t3]; dec(i) != pos[t1]; i = dec(i)) {
            std::swap(pos[t[i]], pos[t[dec(i)]]);
            std::swap(t[i], t[dec(i)]);
        }
    }
};

class TSP {
    time_point<system_clock> end;
    std::mt19937 rng;

    Matrix<weight_t> d;
    Tour tour, best;

    std::vector<std::vector<index_t>> candidates;

    weight_t swap3_dispatch(
        index_t t1, index_t t2,
        index_t t3, index_t t4,
        index_t t5, index_t t6,
        size_t dir  // 0: next, 1: prev
    ) {
        weight_t gain = d(t1, t2) + d(t3, t4) + d(t5, t6)
             - d(t1, t6) - d(t2, t3) - d(t4, t5);
        if (gain <= 0) return 0;
        for (index_t x : {t1, t2, t3, t4, t5, t6}) tour.reset_dlb(x);

        if (dir) {
            tour.reverse(t6, t4);
            tour.reverse(t3, t1);
        } else {
            std::vector<index_t> t;
            t.reserve(n());

            for (size_t i = tour.pos[t2]; i != tour.pos[t5]; i = tour.inc(i)) {
                t.push_back(tour.t[i]);
            }
            t.push_back(t5);
            for (size_t i = tour.pos[t4]; i != tour.pos[t1]; i = tour.inc(i)) {
                t.push_back(tour.t[i]);
            }
            t.push_back(t1);

            size_t i4 = tour.pos[t4];
            for (size_t i = 0; i < t.size(); ++i) {
                tour.t[(i4 + i) % n()] = t[i];
                tour.pos[t[i]] = (i4 + i) % n();
            }
        }
        return gain;
    }

    weight_t opt_2h() {
        weight_t gain = 0;

        bool improved = 1;
        while (improved && time() > 0ns) {
            improved = 0;

            for (index_t t1 : tour.t) { // For some reason this is better?
                if (tour.dlb[t1].try_skip_2()) continue;

                index_t t2 = tour.next(t1);

                for (index_t t3 : candidates[t1]) {
                    weight_t g = d(t1, t2) - d(t1, t3);
                    if (g <= 0) break;  // For faster convergence?

                    index_t t4 = tour.next(t3);
                    if (t1 == t4 || t2 == t3) continue;
                    
                    g += d(t3, t4) - d(t2, t4);
                    if (g > 0) {
                        /*
                        if (tour.segment_length(t2, t3) < n() / 2) {
                            tour.reverse(t2, t3);
                        } else {
                            tour.reverse(t4, t1);
                        }
                        */
                        tour.reverse(t2, t3);
                        for (index_t x : {t1, t2, t3, t4}) tour.reset_dlb(x);

                        gain += g;
                        improved = 1;
                        break;
                    }

                    // variant B
                    index_t b = tour.next(t2);
                    if (b != t3) {
                        g = d(t1, t2) + d(t2, b) + d(t3, t4)
                            - d(t1, b) - d(t3, t2) - d(t2, t4);
                        if (g > 0) {
                            tour.swap2b(t1, t2, t3, t4, b);
                            for (index_t x : {t1, t2, t3, t4, b}) tour.reset_dlb(x);
                            gain += g;
                            improved = 1;
                            break;
                        }
                    }

                    // variant C
                    index_t c = tour.prev(t3);
                    if (c != t2) {
                        g = d(t1, t2) + d(t3, c) + d(t3, t4)
                            - d(t1, t3) - d(t2, t3) - d(t4, c);
                        if (g > 0) {
                            tour.swap2c(t1, t2, t3, t4, c);
                            for (index_t x : {t1, t2, t3, t4, c}) tour.reset_dlb(x);
                            gain += g;
                            improved = 1;
                            break;
                        }
                    }
                }
            }
        }
        return gain;
    }

    weight_t opt_3() {
        weight_t gain = 0;

        bool improved = 1;
        while (improved && time() > 0ns) {
            improved = 0;

            for (index_t t1 : tour.t) {
                if (tour.dlb[t1].try_skip_3()) continue;

                index_t t2 = tour.next(t1);
                size_t i2 = tour.pos[t2];
                for (index_t t3 : candidates[t2]) {
                    size_t i3 = tour.pos[t3];

                    if (d(t1, t2) <= d(t2, t3)) break;
                    if (t1 == t3) continue;
                    
                    for (size_t dir = 0; dir <= 1; ++dir) {
                        index_t t4 = dir ? tour.prev(t3) : tour.next(t3);
                        if (t1 == t4) continue;
                        for (index_t t5 : candidates[t4]) {
                            if (in_list(t5, t1, t2, t3)) continue;

                            size_t i5 = tour.pos[t5];
                            if (!(
                                (i5 < i3 && i3 < i2) ||
                                (i3 < i2 && i2 < i5) ||
                                (i2 < i5 && i5 < i3)
                            )) continue;

                            if (d(t3, t4) <= d(t2, t3) - d(t1, t2) + d(t4, t5)) break;

                            index_t t6 = tour.next(t5);
                            size_t i6 = tour.pos[t6];

                            weight_t g = swap3_dispatch(t1, t2, t3, t4, t5, t6, dir);
                            if (g > 0) {
                                gain += g;
                                improved = 1;
                                goto outer;
                            }
                        }
                    }
                }
outer:
                continue;
            }
        }
        return gain;
    }

public:
    static constexpr size_t LIMIT = 20;

    TSP(const std::vector<Point>& p, nanoseconds duration, uint_fast32_t seed = -88888888):
        end(system_clock::now() + duration - 5ms), rng(seed)
    {
        size_t n = p.size();

        d = Matrix<weight_t>(n, n); 
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                d(i, j) = i == j ? INT_MAX : std::round(p[i].d(p[j]));
            }
        }

        // For equivalence to benchmark, switch to heap
        candidates.resize(n);
        size_t m = std::min(n, LIMIT);
        for (index_t i = 0; i < n; ++i) {
            auto& c = candidates[i];
            c.resize(n);
            std::iota(c.begin(), c.end(), 0);
            std::partial_sort(c.begin(), c.begin() + m, c.end(),
                [&](index_t a, index_t b){ return d(i, a) < d(i, b); });
            c.resize(m);
        }
    }

    size_t n() const {
        return d.n();
    }

    weight_t cost(const Tour& tour) const {
        weight_t c = d(tour.t[n() - 1], tour.t[0]);
        for (size_t i = 0; i < n() - 1; ++i) {
            c += d(tour.t[i], tour.t[i + 1]);
        }
        return c;
    }

    nanoseconds time() const {
        return end - system_clock::now();
    }

    Tour gnn() const {
        std::vector<index_t> v(n());
        v[0] = 0;

        std::vector<bool> used(n(), 0);
        used[0] = 1;

        for (size_t i = 1; i < n(); ++i) {
            index_t best = -1;
            for (size_t j = 0; j < n(); ++j) {
                if (!used[j] && (best == -1 ||
                    d(v[i - 1], j) < d(v[i - 1], best))) {
                    best = j;
                }
            }
            v[i] = best;
            used[best] = 1;
        }

        return Tour(std::move(v));
    }

    Tour solve() {
        tour = gnn();
        best = tour;
        weight_t cost_tour = cost(tour);
        weight_t cost_best = cost_tour;
        
        while (time() > 0ns) {
            tour.perturb(rng);

            // If we recalculate at the end, need not use retvals of opt*
            opt_2h();
            opt_3();
            cost_tour = cost(tour);

            // fprintf(stderr, "%ld %ld\n", cost_tour, cost_best);

            if (cost_tour < cost_best) {
                cost_best = cost_tour;
                best = tour;
            } else {
                // cost_tour = cost_best;
                tour = best;
            }
        }

        // fprintf(stderr, "%ld\n", cost_best);
        return best;
    }
};

#ifdef LOCAL
/*
    For 35/50, need OPT~1.02 on average
    OPT<1.02 for n<=300
    OPT=1.05 for d657.tsp
    OPT=1.08 for dsj1000.tsp (OPT=1.04 if 60s instead)

    Scanner is crudely written. You may have to manually adjust tsp files.
*/
const std::unordered_map<std::string, weight_t> optimal = {
    {"berlin52.tsp", 7542},
    {"d198.tsp", 15780},
    {"d657.tsp", 48912},
    {"dsj1000.tsp", 18660188},
    {"lin105.tsp", 14379},
    {"lin318.tsp", 42029},
    {"eil101.tsp", 629},
    {"kroA200.tsp", 29368},
    {"pr152.tsp", 73682},
    {"sample.tsp", 276},
    {"lu980.tsp", 11340},
    {"qa194.tsp", 9352},
    {"uy734.tsp", 79114},
    {"wi29.tsp", 27603},
    {"zi929.tsp", 95345},
};

char *gnu_basename(char *path) {
    char *base = strrchr(path, '/');
    return base ? base + 1 : path;
}

double tsplib(char* name) {
    FILE* pf = fopen(name, "r");
    char buf[256];

    auto it = optimal.find(basename(name));
    if (it == optimal.end()) return 0;
    double c0 = it->second;

    while (1) {
        fgets(buf, 256, pf);
        if (strncmp(buf, "NODE_COORD_SECTION", 18) == 0) break;
    }
    std::vector<Point> p;
    while (1) {
        Point pt;
        if (fscanf(pf, "%*d %lf %lf", &pt.x, &pt.y) != 2) break;
        p.push_back(pt);
    }
    
    TSP tsp(p, 2s);
    Tour t1 = tsp.gnn();
    Tour t2 = tsp.solve();
    double c1 = tsp.cost(t1);
    double c2 = tsp.cost(t2);
    double score = pow(0.02, (c2 - c0) / (c1 - c0));
    printf("%32s: %9.0lf (%.3lfOPT) -> %9.0lf (%.3lfOPT) (score=%.2lf)\n", name, c1, c1 / c0, c2, c2 / c0, score);
    return score;
}
#endif

// Do not modify
void kattis() {
    size_t n;
    scanf("%ld", &n);
    std::vector<Point> p(n);
    for (auto& [x, y] : p) scanf("%lf %lf", &x, &y);
    TSP tsp(p, 2s);
    Tour tour = tsp.solve();
    for (auto x : tour.vec()) printf("%ld\n", x);
}

int main(int argc, char** argv) {
    if (argc > 1) {
#ifdef LOCAL
        double total = 0;
        double read = argc - 1;
        for (size_t i = 1; i < argc; ++i) {
            double score = tsplib(argv[i]);
            if (score == 0) read -= 1;
            total += score;
        }
        printf("mean: %lf\n", total / read);
#endif
    } else {
        kattis();
    }
}
