# TPGO-RS
Translation-Only Pose Graph Optimization through a Riemannian Staircase

C++ Project

Dependencies:
- [ROFL](https://github.com/dlr1516/rofl) library
- Riemannian Optimization library [ROPTLIB](https://github.com/whuang08/ROPTLIB) by Wen Huang et al. hard-packaged locally
- qr_unique() function from [MANOPT](https://www.manopt.org/) by Boumal et al. compiled through MATLAB coder and hard-packaged locally (in include/thirdparty/qr_unique_sizeless) for increased accuracy in retractions on Stiefel manifold

