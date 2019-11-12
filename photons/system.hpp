#pragma once
#include <cmath>
#include <random>
#include <vector>

namespace photons {

struct system {
  using real_type = float;

  system() = default;
  explicit system(int n) : pos_x(n), pos_y(n), v_x(n), v_y(n), weights(n, 1) {}

  auto size() const noexcept { return pos_x.size(); }

  std::vector<real_type> pos_x{};
  std::vector<real_type> pos_y{};
  std::vector<real_type> v_x{};
  std::vector<real_type> v_y{};
  std::vector<real_type> weights{};
};

template <typename RNG>
inline void generate_random(system& sys, RNG& rng) noexcept {
  constexpr float two_pi = 6.283185307;
  std::uniform_real_distribution<system::real_type> dist{-1, 1};
  for (int i = 0; i < sys.size(); ++i) {
    const auto radius = 0.2 * std::sqrt(dist(rng));
    const auto angle = two_pi * dist(rng);
    sys.pos_x[i] = radius * std::cos(angle);
    sys.pos_y[i] = radius * std::sin(angle);
    // sys.pos_x[i] = 0.5 * dist(rng);
    // sys.pos_y[i] = 0 * dist(rng);
    // sys.v_x[i] = dist(rng);
    // sys.v_y[i] = dist(rng);
    sys.v_x[i] = 0;
    sys.v_y[i] = -1;
    const auto v_inv_norm =
        1.0f / std::sqrt(sys.v_x[i] * sys.v_x[i] + sys.v_y[i] * sys.v_y[i]);
    sys.v_x[i] *= v_inv_norm;
    sys.v_y[i] *= v_inv_norm;
  }
}

template <typename RNG>
inline void advance(system& sys, RNG& rng) noexcept {
  constexpr float two_pi = 6.283185307;
  // constexpr float g = 0.9;
  constexpr auto square = [](auto x) { return x * x; };
  std::uniform_real_distribution<system::real_type> dist{0, 1};

  for (int i = 0; i < sys.size(); ++i) {
    const auto time = -std::log(dist(rng)) * 1e-1;
    sys.pos_x[i] += time * sys.v_x[i];
    sys.pos_y[i] += time * sys.v_y[i];
    const float absorption = 0.01;
    sys.weights[i] *= (1 - absorption);
    if (dist(rng) > absorption) continue;
    // const auto angle = two_pi * dist(rng);
    // const auto cos_angle = std::cos(angle);
    // const auto sin_angle = std::sin(angle);
    // const auto cos_angle =
    //     (1 + square(g) - square((1 - g * g) / (1 - g + 2 * g * dist(rng)))) /
    //     (2 * g);
    const auto cos_angle = 2 * dist(rng) - 1;
    const auto sin_angle =
        (2 * (dist(rng) < 0.5) - 1) * std::sqrt(1 - square(cos_angle));
    const auto new_x = cos_angle * sys.v_x[i] - sin_angle * sys.v_y[i];
    const auto new_y = sin_angle * sys.v_x[i] + cos_angle * sys.v_y[i];
    // const auto new_x =
    //     cos_angle * sys.v_x[i] - sin_angle * std::sqrt(1 -
    //     square(sys.v_x[i]));
    // const auto new_y =
    //     cos_angle * sys.v_y[i] +
    //     sin_angle * sys.v_y[i] * sys.v_x[i] / std::sqrt(1 -
    //     square(sys.v_x[i]));
    sys.v_x[i] = new_x;
    sys.v_y[i] = new_y;
  }
}

namespace experimental {

template <typename RNG>
inline void advance(system& sys, RNG& rng) noexcept {
  constexpr float two_pi = 6.283185307;
  // constexpr float g = 0.9;
  constexpr auto square = [](auto x) { return x * x; };
  std::uniform_real_distribution<system::real_type> dist{0, 1};

  constexpr auto n = [](float x, float y) { return 1 + (1 - 0.5 * x * x); };
  constexpr auto grad_n = [](float x, float y) { return -1.0 * x; };

  for (int i = 0; i < sys.size(); ++i) {
    // const auto time = -std::log(dist(rng)) * 1e-2;
    const auto time = 1e-2f;
    const auto inv_n = 1 / n(sys.pos_x[i], sys.pos_y[i]);
    const auto grad = grad_n(sys.pos_x[i], sys.pos_y[i]);
    sys.pos_x[i] += time * sys.v_x[i] * inv_n;
    sys.pos_y[i] += time * sys.v_y[i] * inv_n;

    sys.v_x[i] += 0.5 * time * grad * inv_n;
    const auto v_inv_norm =
        1.0f / std::sqrt(sys.v_x[i] * sys.v_x[i] + sys.v_y[i] * sys.v_y[i]);
    sys.v_x[i] *= v_inv_norm;
    sys.v_y[i] *= v_inv_norm;

    const float absorption = 0.0001;
    const float scattering = 0.95;
    if (dist(rng) > absorption) continue;
    if (dist(rng) > scattering) {
      sys.weights[i] = 0;
      continue;
    }
    const auto angle = two_pi * dist(rng);
    const auto cos_angle = std::cos(angle);
    const auto sin_angle = std::sin(angle);
    // const auto cos_angle =
    //     (1 + square(g) - square((1 - g * g) / (1 - g + 2 * g * dist(rng)))) /
    //     (2 * g);
    // const auto cos_angle = 2 * dist(rng) - 1;
    // const auto sin_angle =
    //     (2 * (dist(rng) < 0.5) - 1) * std::sqrt(1 - square(cos_angle));
    const auto new_x = cos_angle * sys.v_x[i] - sin_angle * sys.v_y[i];
    const auto new_y = sin_angle * sys.v_x[i] + cos_angle * sys.v_y[i];
    // const auto new_x =
    //     cos_angle * sys.v_x[i] - sin_angle * std::sqrt(1 -
    //     square(sys.v_x[i]));
    // const auto new_y =
    //     cos_angle * sys.v_y[i] +
    //     sin_angle * sys.v_y[i] * sys.v_x[i] / std::sqrt(1 -
    //     square(sys.v_x[i]));
    sys.v_x[i] = new_x;
    sys.v_y[i] = new_y;
  }
}

}  // namespace experimental

}  // namespace photons