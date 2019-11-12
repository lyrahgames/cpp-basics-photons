#include <cmath>
#include <photons/application.hpp>

using namespace std;

namespace photons {

application::application() : window{{500, 500}, "Photons"}, sys{1000} {
  view.setCenter(0, 0);

  photons::generate_random(sys, rng);
}

void application::resize(int w, int h) {
  view.setSize(w, h);
  view.zoom(fov / h);
  window.setView(view);
}

void application::render() {
  const float radius = fov / 500;
  sf::CircleShape circle{radius};
  circle.setOrigin(radius, radius);

  for (int i = 0; i < sys.size(); ++i) {
    circle.setPosition(sys.pos_x[i], sys.pos_y[i]);
    circle.setFillColor(sf::Color{255 * sys.weights[i], 255 * sys.weights[i],
                                  255 * sys.weights[i], 128});
    window.draw(circle);
  }
}

void application::execute() {
  bool update = false;
  sf::Vector2i old_mouse{};

  while (window.isOpen()) {
    const auto mouse = sf::Mouse::getPosition(window);

    sf::Event event{};
    while (window.pollEvent(event)) {
      switch (event.type) {
        case sf::Event::Closed:
          window.close();
          break;

        case sf::Event::Resized:
          resize(event.size.width, event.size.height);
          break;

        case sf::Event::MouseWheelMoved:
          fov *= exp(-event.mouseWheel.delta * 0.05f);
          fov = clamp(fov, 1e-6f, 1000.f);
          update = true;
          break;

        case sf::Event::KeyPressed:
          switch (event.key.code) {
            case sf::Keyboard::Escape:
              window.close();
              break;
          }
          break;
      }
    }

    if (sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
      const auto mouse_move =
          window.mapPixelToCoords(mouse) - window.mapPixelToCoords(old_mouse);
      view.move(-mouse_move);
      update = true;
    }

    if (update) {
      resize(window.getSize().x, window.getSize().y);
      update = false;
    }

    photons::experimental::advance(sys, rng);

    window.clear();
    render();
    window.display();

    old_mouse = mouse;
  }
}

}  // namespace photons