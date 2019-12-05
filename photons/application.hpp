#pragma once
#include <SFML/Graphics.hpp>
#include <photons/system.hpp>
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
  std::mt19937 rng{std::random_device{}()};
};

}  // namespace photons