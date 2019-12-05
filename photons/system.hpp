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
    sys.pos_x[i] = 10 * radius * std::cos(angle);
    sys.pos_y[i] = 0.5 * radius * std::sin(angle);
    float cos_angle = cos(-two_pi / 8);
    float sin_angle = sin(-two_pi / 8);
    const auto new_x = cos_angle * sys.pos_x[i] - sin_angle * sys.pos_y[i];
    const auto new_y = sin_angle * sys.pos_x[i] + cos_angle * sys.pos_y[i];
    sys.pos_x[i] = new_x - 10.0f;
    sys.pos_y[i] = new_y + 10.0f;
    // sys.pos_x[i] = 0 + 0.1 * dist(rng);
    // sys.pos_y[i] = 2 + 0.1 * dist(rng);
    // sys.pos_x[i] = 0.5 * dist(rng);
    // sys.pos_y[i] = 0 * dist(rng);
    // sys.v_x[i] = dist(rng);
    // sys.v_y[i] = dist(rng);
    sys.v_x[i] = 1;
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

namespace phase_function {

template <typename RNG>
inline void advance(system& sys, RNG& rng) noexcept {
  constexpr float two_pi = 6.283185307;
  constexpr auto square = [](auto x) { return x * x; };
  std::uniform_real_distribution<system::real_type> dist{0, 1};

  for (int i = 0; i < sys.size(); ++i) {
    const auto time = 1e-3f;
    // const float g = std::clamp(sys.pos_y[i], -1.0f, 1.0f);
    const float g = 0.5f;
    // const float g =
    //     (square(sys.pos_x[i]) + square(sys.pos_y[i]) < 1) ? (-0.8) : (1);
    sys.pos_x[i] += time * sys.v_x[i];
    sys.pos_y[i] += time * sys.v_y[i];

    const float absorption = 0.5;
    const float scattering = 0.95;
    if (dist(rng) < std::exp(-absorption * time)) continue;
    if (dist(rng) > scattering) {
      sys.weights[i] = 0;
      continue;
    }
    const auto cos_angle =
        (1 + square(g) - square((1 - g * g) / (1 - g + 2 * g * dist(rng)))) /
        (2 * g);
    const auto sin_angle =
        ((dist(rng) > 0.5) ? 1 : -1) * std::sqrt(1 - square(cos_angle));
    const auto new_x = cos_angle * sys.v_x[i] - sin_angle * sys.v_y[i];
    const auto new_y = sin_angle * sys.v_x[i] + cos_angle * sys.v_y[i];
    sys.v_x[i] = new_x;
    sys.v_y[i] = new_y;
  }
}

}  // namespace phase_function

namespace optics {

template <typename RNG>
inline void advance(system& sys, RNG& rng) noexcept {
  constexpr auto time = 0.5e-1f;
  std::uniform_real_distribution<system::real_type> dist{0, 1};
  constexpr auto square = [](auto x) { return x * x; };

  const float plane_n_x = 0;
  const float plane_n_y = -1;
  const float plane_d = 0;
  const float plane_t = 1.5;
  const float plane_t2 = plane_t * plane_t;

  for (int i = 0; i < sys.size(); ++i) {
    // t for plane
    const float dot_n_p = plane_n_x * sys.pos_x[i] + plane_n_y * sys.pos_y[i];
    float dot_n_v = plane_n_x * sys.v_x[i] + plane_n_y * sys.v_y[i];
    const float t = (-plane_d - dot_n_p) / dot_n_v;
    if ((0 <= t) && (t <= time)) {
      const float reflexion = 0.5;
      if (dist(rng) < reflexion) {
        // reflexion
        float new_v_x = sys.v_x[i] - 2 * dot_n_v * plane_n_x;
        float new_v_y = sys.v_y[i] - 2 * dot_n_v * plane_n_y;
        sys.pos_x[i] += t * sys.v_x[i];
        sys.pos_y[i] += t * sys.v_y[i];
        sys.v_x[i] = new_v_x;
        sys.v_y[i] = new_v_y;
        sys.pos_x[i] += (time - t) * sys.v_x[i];
        sys.pos_y[i] += (time - t) * sys.v_y[i];
      } else {
        // refraction
        float eta = plane_t;
        if (dot_n_v > 0) {
          eta = -1 / eta;
          // dot_n_v = -dot_n_v;
        }
        const float eta2 = eta * eta;
        const float sqnorm_v =
            sys.v_x[i] * sys.v_x[i] + sys.v_y[i] * sys.v_y[i];
        const float control = sqnorm_v - eta2 * (sqnorm_v - dot_n_v * dot_n_v);
        if (control <= 0) {
          float new_v_x = sys.v_x[i] - 2 * dot_n_v * plane_n_x;
          float new_v_y = sys.v_y[i] - 2 * dot_n_v * plane_n_y;
          sys.v_x[i] = new_v_x;
          sys.v_y[i] = new_v_y;
        } else {
          const float tmp = eta * std::sqrt(control);
          const float new_v_x =
              eta2 * (sys.v_x[i] - dot_n_v * plane_n_x) - tmp * plane_n_x;
          const float new_v_y =
              eta2 * (sys.v_y[i] - dot_n_v * plane_n_y) - tmp * plane_n_y;
          sys.pos_x[i] += t * sys.v_x[i];
          sys.pos_y[i] += t * sys.v_y[i];
          sys.v_x[i] = new_v_x;
          sys.v_y[i] = new_v_y;
          sys.pos_x[i] += (time - t) * sys.v_x[i];
          sys.pos_y[i] += (time - t) * sys.v_y[i];
        }
      }
    } else {
      sys.pos_x[i] += time * sys.v_x[i];
      sys.pos_y[i] += time * sys.v_y[i];

      bool inside = plane_d + dot_n_p > 0;
      const float g = (inside) ? (-0.1) : (0.999);
      const float absorption = (inside) ? (0.01) : (0);
      const float scattering = (inside) ? (1) : (1);
      if (dist(rng) < std::exp(-absorption * time)) continue;
      if (dist(rng) > scattering) {
        sys.weights[i] = 0;
        continue;
      }
      const auto cos_angle =
          (1 + square(g) - square((1 - g * g) / (1 - g + 2 * g * dist(rng)))) /
          (2 * g);
      const auto sin_angle =
          ((dist(rng) > 0.5) ? 1 : -1) * std::sqrt(1 - square(cos_angle));
      const auto new_x = cos_angle * sys.v_x[i] - sin_angle * sys.v_y[i];
      const auto new_y = sin_angle * sys.v_x[i] + cos_angle * sys.v_y[i];
      sys.v_x[i] = new_x;
      sys.v_y[i] = new_y;
    }
  }
}

}  // namespace optics

}  // namespace photons