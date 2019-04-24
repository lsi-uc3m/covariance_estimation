# DCE: Online Drift Covariance Estimation

### Package Summary ###

This repository is for estimating the covariance of drift suffering proprioceptive sensors, using an exteroceptive sensor with known uncertainty. This package is a generic package, which is capable of estimating the covariance of any proprioceptive sensor (any dimension).

**Maintainer Status**: maintained

**Maintainer**: Mostafa Osman <mostafaosman144 AT gmail DOT com>

**Author**: Mostafa Osman, Ahmed Hussein and Abdulla Al-Kaff
  
### Subscribed Topics ###

Proprioceptive sensor odometry with unknown covariance topic

```
/odom_with_no_covariance_topic (nav_msgs/Odometry)
```

Exteroceptive sensor odometry with known covariance topic

```
/odom_with_covariance_topic (nav_msgs/Odometry)
```

### Published Topics ###

Proprioceptive sensor odometry with estimated covariance topic

```
/odom_with_calculated_cov_topic (nav_msgs/Odometry)
```

### Publications ###

M. Osman, A. Hussein, A. Al-Kaff, F. Garcia, and J.M. Armingol (2018). “**Online adaptive covariance estimation approach for multiple odometry sensors fusion**”. *In proceedings of Intelligent Vehicles Symposium (IV2018)*, pp. 1–6. IEEE.
