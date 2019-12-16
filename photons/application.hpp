#pragma once
#include <SFML/Graphics.hpp>
#include <photons/system.hpp>
#include <pxart/pxart.hpp>
#include <random>

namespace photons {

class application {
 public:
  application();

  void execute();

 private:
  void resize(int w, int h);
  void render();

 private:
  sf::RenderWindow window;
  sf::View view;

  float fov = 100;

  system sys;
  system tmp_sys;
  // std::mt19937 rng{std::random_device{}()};
  pxart::mt19937 rng{std::random_device{}};
  // pxart::simd256::mt19937 vrng{std::random_device{}};
  pxart::simd256::xrsr128p vrng{std::random_device{}};
};

}  // namespace photons