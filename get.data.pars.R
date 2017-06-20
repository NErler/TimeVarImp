#########################################
# define parameters for data simulation #
#########################################

# linear sequential ------------------------------------------------------------
if(sim.data.type %in% c("sim.lin.seq")){
  beta <- list(beta.y = c(Intercept = -1, x = 1.2, time = 0.8),
               beta.x = c(Intercept = 2, time = -0.3))
  D <- list(D.y = matrix(c(1.413,  0.166,  0.166, 0.106), 2, 2),
            D.x = matrix(c(1.500, -0.202, -0.202, 0.165), 2, 2))

  beta.analysis <- beta$beta.y
  Sigma = NULL
  sig <- c(sig.y = 0.4, sig.x = 0.3)
}

if(sim.data.type %in% c("sim.lin.seq.b")){
  beta <- list(beta.y = c(Intercept = -1, x = 1.2, time = 0.8),
               beta.x = c(Intercept = 2, time = -0.3))
  
  D <- matrix(c( 1.413,  0.166,  0.022, -0.015,
                 0.166,  0.106,  0.166, -0.083,
                 0.022,  0.166,  1.500, -0.202,
                 -0.015, -0.083, -0.202,  0.165), 4, 4)
  
  beta.analysis <- beta$beta.y
  Sigma = NULL
  sig <- c(sig.y = 0.4, sig.x = 0.3)
}


# quadratic sequential ---------------------------------------------------------
if(sim.data.type %in% c("sim.qdr.seq")){
  beta <- list(beta.y = c(Intercept = -1, x = 1.2, x2 = -0.35, time = 0.8),
               beta.x = c(Intercept = 2, time = -0.3))

  D <- list(D.y = matrix(c(1.413,  0.166,  0.166, 0.106), 2, 2),
            D.x = matrix(c(1.500, -0.202, -0.202, 0.165), 2, 2))
  beta.analysis <- beta$beta.y
  sig <- c(sig.y = 0.4, sig.x = 0.3)
  Sigma = NULL
}

if(sim.data.type %in% c("sim.qdr.seq.b")){
  beta <- list(beta.y = c(Intercept = -1, x = 1.2, x2 = -0.35, time = 0.8),
               beta.x = c(Intercept = 2, time = -0.3))

  D <- matrix(c( 1.413,  0.166,  0.022, -0.015,
                 0.166,  0.106,  0.166, -0.083,
                 0.022,  0.166,  1.500, -0.202,
                -0.015, -0.083, -0.202,  0.165), 4, 4)

  beta.analysis <- beta$beta.y
  sig <- c(sig.y = 0.4, sig.x = 0.3)
  Sigma = NULL
}


# multivariate normal ----------------------------------------------------------
if(sim.data.type %in% c("sim.lin.joint.b")){
  beta <- list(beta.y = c(Intercept = -1, time = 0.8),
               beta.x = c(Intercept = 2, time = -0.3))

    D <- matrix(c( 1.413,  0.166,  0.022, -0.015,
                 0.166,  0.106,  0.166, -0.083,
                 0.022,  0.166,  1.500, -0.202,
                -0.015, -0.083, -0.202,  0.165), 4, 4)

  Sigma = NULL
  sig <- c(sig.y = 0.4, sig.x = 0.3)
}


if(sim.data.type %in% c("sim.lin.joint.be")){
  beta <- list(beta.y = c(Intercept = -1, time = 0.8),
               beta.x = c(Intercept = 2, time = -0.3))
  D <- matrix(c( 1.413,  0.166,  0.022, -0.015,
                 0.166,  0.106,  0.166, -0.083,
                 0.022,  0.166,  1.500, -0.202,
                -0.015, -0.083, -0.202,  0.165), 4, 4)

  Sigma <- matrix(c(0.4, 0.2, 0.2, 0.3), 2,2)
  sig <- c(sig.y = 0.4, sig.x = 0.3)
}


# ------------------------------------------------------------------------------
const <- list("sim.lin.seq"        = c(-3.2, -0.5),
              "sim.lin.seq.b"      = c(-3.2, -0.5),
              "sim.qdr.seq"        = c(-1.85, -0.7),
              "sim.qdr.seq.b"      = c(-1.85, -0.7),
              "sim.lin.joint.b"    = c(-1.6,  -0.5),
              "sim.lin.joint.be"   = c(-1.6,  -0.5))

cutoffs <- c("sim.lin.seq"       =  0.0,
             "sim.lin.seq.b"     =  0.0,
             "sim.qdr.seq"       = -0.3,
             "sim.qdr.seq.b"     = -0.3,
             "sim.lin.joint.b"   =  -0.3,
             "sim.lin.joint.be"  =  -0.3)
