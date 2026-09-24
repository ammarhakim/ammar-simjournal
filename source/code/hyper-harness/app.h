#include <mug.h>

//
// Hyperbolic solver app
//

// typedefs for various functions needed in eqn_sys: one pointer for
// each direction
typedef void (*eqn_flux_t)(void *ctx, const double *q, double *fout);
typedef void (*eqn_projon_left_ev_t)(void *ctx, const double *q, const double *vin, double *vout);
typedef void (*eqn_recwith_right_ev_t)(void *ctx, const double *q, const double *vin, double *vout);
typedef void  (*eqn_rescale_right_ev_t)(void *ctx, const double *q, const double *w, double *revout);
typedef double (*eqn_ev_t)(void *ctx, const double *q, double *ev);
typedef double (*fluct_t)(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq);
typedef void (*rotate_to_local_t)(void *ctx, double n[3], double tau1[3], double tau2[3], const double *qloc, double *qglo);
typedef void (*rotate_to_global_t)(void *ctx, double n[3], double tau1[3], double tau2[3], const double *qglo, double *qloc);

struct eqn_sys {
  int meqn, mtotal;
  void *eqn_ctx;

  eqn_flux_t flux[GKYL_MAX_DIM];
  eqn_projon_left_ev_t projon_left_ev[GKYL_MAX_DIM];
  eqn_recwith_right_ev_t recwith_right_ev[GKYL_MAX_DIM];
  eqn_rescale_right_ev_t rescale_right_ev[GKYL_MAX_DIM];
  eqn_ev_t ev[GKYL_MAX_DIM];
  fluct_t fluct[GKYL_MAX_DIM];
  rotate_to_local_t rotate_to_local;
  rotate_to_global_t rotate_to_global;
};

typedef struct hyper_app hyper_app;

enum hyper_scheme_type { SCHEME_FD, SCHEME_FV };

// Base reconstruction scheme to use
enum mp_recon {
  MP_U3, // upwind-biased 3rd order
  MP_U5, // upwind-biased 5th order
  MP_U1, // upwind-biased 1st order
  MP_C2, // centered second-order
  MP_C4, // centered fourth-order
  MP_C6, // centered sixth-order
};

struct hyper_app_inp {
  const char name[128];

  int nframe;
  double tend;

  int ndim;
  int cells[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  enum hyper_scheme_type scheme_type;
  struct eqn_sys hyper_eqn;
  bool has_diffusive_eqn;
  struct eqn_sys diff_eqn;

  enum mp_recon recon_type;
  bool use_char_limiters;
  
  evalf_t init;
  void *init_ctx;
};

hyper_app *hyper_app_new(struct hyper_app_inp *inp);
struct gkyl_update_status hyper_app_update(hyper_app *app, double dt);
void hyper_app_release(hyper_app *app);
