#include <gkylzero.h>

typedef struct adiff_app adiff_app;

enum adiff_advection_scheme {
  SCHEME_C2_C2, // second-order central for both advection and diffusion
  SCHEME_U1_C2, // first-order upwind for advection, second-order central for diffusion
  SCHEME_U3_C4, // second-order upwind for advection, fourth-order central for diffusion
};  

struct adiff_app_inp {
  char name[128]; // name of simulation
  
  int ndim; // number of dimensions (defaults to 2)
  int cells[3]; // cells in each direction
  double lower[3]; // lower extents of domain
  double upper[3]; // upper extents of domain

  int nframe; // number of output frames to write
  double tend; // time to run simulation
  double cfl_frac; // Adjust CFL number by this (< 1.0)

  enum adiff_advection_scheme scheme; // scheme to use for advection term

  void *init_ctx; // context object for ICs
  evalf_t init; // initial conditions

  void *velocity_ctx; // context object for advection velocity
  evalf_t velocity; // function for advection velocity

  double alpha; // diffusion coefficient (constant)

  // for debugging only  
  int max_steps;
};

/**
 * Construct a new advection-diffusion app object.
 *
 * @param inp Input to construct app
 * @return new advection-diffusion app
 */
adiff_app *adiff_app_new(const struct adiff_app_inp *inp);

/**
 * Run simulation
 *
 * @param app App object
 */
void adiff_app_run(adiff_app *app);

/**
 * Release memory for advection-diffusion app.
 *
 * @param app App object to release
 */
void adiff_app_release(adiff_app *app);
