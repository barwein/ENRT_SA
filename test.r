
source("src/sensitivity_params.r")
source("src/bias_adjustment.r")

n_e <- 100
n_a <- 200
pz <- 0.5

X_e <- matrix(rbinom(n=n_e*3, size = 1, prob = 0.4), nrow=n_e, ncol = 3)
X_a <- matrix(rbinom(n=n_a*3, size = 1, prob = 0.3), nrow = n_a, ncol = 3)

ego_index <- rep(seq(n_e), each = 2)

rho_vec <- seq(1e-3,1e-2, 1e-3)
m_vec_ee <- seq(1,20,2)
m_vec_ae <- seq(1,20,2)

# EGO-EGO
pi_homo_ee <- pi_homo(rho_vec = rho_vec,
                      n_e = n_e,
                      pz = pz,
                      type = "ego")

pi_num_homo_ee <- pi_homo(m_vec = m_vec_ee,
                          n_e = n_e,
                          pz = pz,
                          type = "ego")

pi_hetero_ee <- pi_hetero(X_e = X_e,
                          gamma = seq(0.1,3,0.25), 
                          dist = "norm", 
                          p = 1,
                          pz = pz)

pi_num_hetero_ee <- pi_hetero(X_e = X_e,
                          gamma = -1,
                          m_vec = m_vec_ee,
                          dist = "norm", 
                          p = 1,
                          pz = pz)


# ALTER-EGO
pi_homo_ae <- pi_homo(rho_vec = rho_vec,
                      n_e = n_e,
                      n_a = n_a,
                      type = "alter",
                      pz = pz)

pi_num_homo_ae <- pi_homo(m_vec = m_vec_ae,
                          n_e = n_e,
                          n_a =n_a,
                          pz = pz,
                          type = "alter")

pi_hetero_ae <- pi_hetero(X_e = X_e,
                          X_a = X_a,
                          gamma = seq(0.1,3,0.25), 
                          dist = "norm",
                          p=1,
                          ego_index = ego_index,
                          pz = pz)

pi_num_hetero_ae <- pi_hetero(X_e=X_e,
                              X_a=X_a,
                              m_vec = m_vec_ae, 
                              gamma=-1,
                              dist = "norm",
                              p=2,
                              ego_index = ego_index,
                              pz=pz)

# TEST BIAS ADJUSTMENT FOR IE
ui_alters <- runif(n_a)
mu_01 <- 1 / (1 + exp(-(ui_alters + 0.5)))
mu_00 <- 1 / (1 + exp(-(ui_alters)))

esti_mat <- as.matrix(cbind(mu_01, mu_00))

ie_rd_naive <- mean(esti_mat %*% c(1,-1))
ie_rr_naive <- sum(esti_mat %*% c(1,0)) / sum(esti_mat %*% c(0,1))

ie_homo_num <- ie_pi_homo_point_grid_(mu_01 =  mu_01,
                                      mu_00 = mu_00,
                                      pi_list = pi_num_homo_ae,
                                      pz = pz)

ie_hetero_num <- ie_pi_homo_point_grid_(mu_01 =  mu_01,
                                        mu_00 = mu_00,
                                        pi_list = pi_num_hetero_ae,
                                        pz = pz)



# TEST BIAS ADJUSTMENT FOR DE
ui_egos <- runif(n_e)
mu_10 <- 1 / (1 + exp(-(ui_egos + 0.5)))
mu_00 <- 1 / (1 + exp(-(ui_egos + 0.1)))
mu_01 <- 1 / (1 + exp(-(ui_egos + 0.25)))

kappa_vec <- seq(1, 1.25, 0.025)

de_hetero_num <- de_grid_multi_pi_kappa(mu_10 = mu_10,
                                        mu_00 = mu_00,
                                        mu_01 = mu_01,
                                        pi_list = pi_num_hetero_ee, 
                                        kappa_vec = kappa_vec)

de_homo_num <- de_grid_multi_pi_kappa(mu_10 = mu_10,
                                        mu_00 = mu_00,
                                        mu_01 = mu_01,
                                        pi_list = pi_num_homo_ee, 
                                        kappa_vec = kappa_vec)


