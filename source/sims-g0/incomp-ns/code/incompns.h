#include <gkylzero.h>

typedef struct incompns_app incompns_app;

// Input parameters to app
struct incompns_inp {
  char name[128]; // name of app: used as output prefix
  
  int ndim; // dimension
  double lower[3], upper[3]; // lower, upper bounds
  int cells[3]; // config-space cells

  double cfl_frac; // CFL fraction to use

  double nu; // dynamic viscosity (units of L^2/T)
  double kappa; // density diffusion (units of L^2/T)

  int num_periodic_dir; // number of periodic directions
  int periodic_dirs[3]; // list of periodic directions

  // non-periodic BCs
  enum gkyl_species_bc_type bcx[2], bcy[2], bcz[2];
  // for driven simulations, set velocity of fluid at boundary
  double bcx_vel[2], bcy_vel[2], bcz_vel[2];

  evalf_t rho_init, vel_init; // density and velocity ICs
  void *ctx; // eval context for initial conditions
};

// Input to the App run method
struct incomp_run_params {
  double tstart, tend; // start and end time of simulation
  int nframe;          // number of frames to write
  bool use_single_step; // should we single-step through the simulation?
  int max_step; // maximum number of steps (only useful in use_single_step = true mode)
};

/**
 * Construct new incompressible NS App
 *
 * @param inp App inputs
 */
incompns_app *incompns_app_new(struct incompns_inp *inp);

/**
 * Apply initial condition.
 *
 * @param app App object
 */
void incompns_app_apply_ic(incompns_app *app);

/**
 * Write data to gkyl files.
 *
 * @param app App object
 * @param tm Time at which data is written
 * @param frame Output frame number
 */
void incompns_write(incompns_app *app, double tm, int frame);

/**
 * Advance simulation by a suggested time-step 'dt'. The dt may be too
 * large in which case method will attempt to take a smaller time-step
 * and also return it as the 'dt_actual' field of the status
 * object. If the suggested time-step 'dt' is smaller than the largest
 * stable time-step the method will use the smaller value instead,
 * returning the larger time-step in the 'dt_suggested' field of the
 * status object. If the method fails to find any stable time-step
 * then the 'success' flag will be set to 0. At that point the calling
 * code must abort the simulation as this signals a catastrophic
 * failure and the simulation can't be safely continued.
 * 
 * @param app App object.
 * @param dt Suggested time-step to advance simulation
 * @return Status of update.
 */
struct gkyl_update_status incompns_update(incompns_app *app, double dt);

/**
 * Run the simulation from @a tstart to @a tend. Returns true if
 * simulation succeeded.
 *
 * @param app App object
 * @param run_param Run parameters
 */
bool incompns_run(incompns_app *app, struct incomp_run_params run_param);

/**
 * Release App object
 *
 * @param app App object to release
 */
void incompns_app_release(incompns_app *app);
