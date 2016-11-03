void cons_qty_init(struct cell_var * cv, const struct flu_var FV);
int cons2prim(struct i_f_var * ifv);
int cons_qty_update
(struct cell_var * cv, const struct mesh_var mv,
 const struct flu_var FV, const double tau);
int cons_qty_update_corr_ave_P
(struct cell_var * cv, const struct mesh_var mv,
 const struct flu_var FV, const double tau);
int cons_qty_update_corr_ave_P_all_wave
(struct cell_var * cv, const struct mesh_var mv,
 const struct flu_var FV, const double tau);

struct cell_var cell_mem_init(const struct mesh_var mv, struct flu_var * FV);
void vol_comp(struct cell_var * cv, const struct mesh_var mv);
void cell_pt_clockwise(struct mesh_var * mv);
void cell_rel(struct cell_var * cv, const struct mesh_var mv);
void cell_centroid(struct cell_var * cv, const struct mesh_var mv);


void slope_limiter(struct cell_var * cv, const struct mesh_var mv, const struct flu_var FV);


void cons_qty_copy_cv2ifv(struct i_f_var *ifv, struct cell_var cv, int c);
void prim_var_copy_ifv2FV(struct i_f_var ifv, struct flu_var *FV, int c);
void flux_copy_ifv2cv(struct i_f_var ifv, struct cell_var *cv, int k, int j);


int fluid_var_update(struct flu_var *FV, struct cell_var cv);
int interface_var_init
(const struct cell_var cv, const struct mesh_var mv, struct i_f_var * ifv,
 struct i_f_var * ifv_R, const int k, const int j, const int i);
double tau_calc(const struct cell_var cv, const struct mesh_var mv);


void Roe_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R);
void HLL_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R);
void Riemann_exact_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R);
void GRP_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R, double tau);


void finite_volume_scheme(struct flu_var *FV, const struct mesh_var mv, const char *scheme);
