#include <gkylzero.h>

typedef struct adiff_app adiff_app;

struct adiff_app_inp {
  int cells[2]; // cells in each direction
  double lower[2]; // lower extents of domain
  double upper[2]; // upper extents of domain

  int nframe; // number of output frames to write
  double tend; // time to run simulation
  double cfl; // CFL condition

  evalf_t init; // initial conditions
  evalf_t vel;  // function for advection velocity
  double alpha; // diffusion coefficient
};

/**
 * Construct a new advection-diffusion app object.
 *
 * @param inp Input to construct app
 * @return new advection-diffusion app
 */
adiff_app* adiff_app_new(const struct adiff_app_inp *inp);

/**
 * Release memory for advection-diffusion app.
 *
 * @param app App object to release
 */
void adiff_app_release(adiff_app *app);
