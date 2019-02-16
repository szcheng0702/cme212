#include "SFML/Graphics.hpp"
#include <cstdlib>

/**
 * @brief Implements utilities to display grid of square cells
 * on screen.
 *
 * Built atop of SFML graphics and animation library. Cells are
 * stored in row-major order in vector container viewgrid_.
 *
 * @tparam I Integer type for indexing cell position in the grid.
 * @tparam T Real floating point type
 */
template <typename I, typename T>
class Viewer
{
private:
  I delta_;      ///< Size of each square cell in pixels
  I xsize_;      ///< Number of cell rows
  I ysize_;      ///< Number of cell columns
  I gridsize_;   ///< Number of cells in grid = xsize_ * ysize_
  std::vector<sf::RectangleShape> viewgrid_; ///< Cells to plot
  sf::RenderWindow* window_;                 ///< Pointer to the veiwer window
  
public:
  /**
   * @brief Viewer constructor
   *
   * Creates an instance of a class that visualizes a rectngular grid
   * of square cells. Creates an SFML window to display cells.
   *
   * @param[in] xsize Number of cell rows
   * @param[in] ysize Number of cell columns
   * @param[in] title Window title
   * @param[in] delta Size of the square cell (length of the edge) in pixels
   */
  Viewer(I xsize=100, I ysize=100, std::string title = "", I delta=4)
  : delta_(delta), 
    xsize_(xsize),
    ysize_(ysize),
    gridsize_(xsize_*ysize_),
    viewgrid_(gridsize_)
  {
    window_ = new sf::RenderWindow(sf::VideoMode(xsize_*delta_, ysize_*delta_), title);
    window_->clear();
    initGrid();
  }

  /**
   * @brief Viewer destructor
   *
   * Deletes SFML window.
   */
  ~Viewer()
  {
    delete window_;
  }

  /**
   * @brief Checks if SFML window is open.
   *
   * @return true if the window is open
   */
  bool isOpen()
  {
    return window_->isOpen();
  }

  /**
   * @brief Updates color of cell at position (i,j) in display grid
   *
   * @param[in] i row index of the cell
   * @param[in] j column index of the cell
   * @param[in] r red color intensity of the cell at (i,j)
   * @param[in] g green color intensity of the cell at (i,j)
   * @param[in] b blue color intensity of the cell at (i,j)
   *
   * @pre _r_, _g_, and _b_ must be real numbers in [0,1)
   * @pre static_cast to sf::Uint8 must be implemented for _r_, _g_, and _b_
   */  
  void updateCell(I i, I j, T r, T g, T b)
  {
    sf::Uint8 red   = static_cast<sf::Uint8>(r*255.0 + 0.5);
    sf::Uint8 green = static_cast<sf::Uint8>(g*255.0 + 0.5);
    sf::Uint8 blue  = static_cast<sf::Uint8>(b*255.0 + 0.5);
    viewgrid_[i*ysize_ + j].setFillColor(sf::Color(red, green, blue));
  }

  /**
   * @brief Initializes grid of cells
   * 
   * @pre viewgrid_ is a vector of sf::RectangleShape objects
   * @pre viewgrid_ has size = (xsize_ * ysize_)
   */
  void initGrid()
  {
    for (I i = 0; i < xsize_; ++i)
      for (I j = 0; j < ysize_; ++j)
	{
	  I n = i*ysize_ + j;
	  float delta = static_cast<float>(delta_);
	  float x = static_cast<float>(i * delta_);
	  float y = static_cast<float>(j * delta_);
	  viewgrid_[n].setSize(sf::Vector2f(delta, delta));
	  viewgrid_[n].setPosition(sf::Vector2f(x, y));
	  viewgrid_[n].setFillColor(sf::Color(0,0,0));
	}
  }
  
  
  /**
   * @brief Redraws the grid of cells on display
   * 
   */
  void refresh()
  {
    for(I n = 0; n < gridsize_; ++n)
    {
      window_->draw(viewgrid_[n]);
    }
    window_->display();
  }
  
  /**
   * @brief Holds on display window waiting for an event
   * 
   */
  void event_loop()
  {
    sf::Event event;
    while (window_->pollEvent(event)) {
      handle_event(event);
    }
    window_->display();
  }

  /**
   * @brief Event handler
   * 
   * This will close window if you click on close window button.
   * Everything else is ignored.
   *
   * @param[in] event SFML event object
   */
  void handle_event(sf::Event event)
  {
    switch (event.type)
    {
      case sf::Event::MouseButtonPressed: 
        break;
      case sf::Event::MouseButtonReleased:
        break;
      case sf::Event::MouseMoved: 
        break;
      case sf::Event::MouseWheelScrolled: 
        break;
      case sf::Event::KeyPressed: 
        break;
      case sf::Event::Resized: 
        break;

      case sf::Event::Closed:
        window_->close();
        break;
      default:
        break;
    }
  }
};
