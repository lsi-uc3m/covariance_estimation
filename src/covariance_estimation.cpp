    #include <ros/ros.h>
  #include "nav_msgs/Odometry.h"
  #include <math.h>
  #include <string>
  #include <eigen3/Eigen/Dense>
  #include <eigen3/unsupported/Eigen/MatrixFunctions>
  #include "XmlRpcValue.h"
  #include "vector"
  #include <iostream>
  #include "tf/transform_datatypes.h"
  #include "tf/transform_broadcaster.h"

  #define STATES_SIZE 6

  ros::Publisher odomWithCalculatedCovPublisher;

  Eigen::VectorXd ExteroceptiveSensorMeasurement_;               //Z
  Eigen::MatrixXd ExteroceptiveSensorCovariance_;                //R

  Eigen::VectorXd ProprioceptiveSensorMeasurement_;              //X
  Eigen::MatrixXd ProprioceptiveSensorCovariance_;               //P
  ros::Time readingTime;

  Eigen::VectorXd StatesConfig_;  //M
  Eigen::VectorXd scaleFactors_;

  int sizeOfTimeBuffer; //size of the innovation memory matrix
  double ExteroCovConfidenceFactor_;
  double MahalanobisThreshold_;
  double C;
  bool firstIteration = true;
  std::string vehicle;
  std::string local_frame_id;

  Eigen::VectorXd loadStatesConfiguration(ros::NodeHandle nh)
  {
    XmlRpc::XmlRpcValue statesConfigXml;
    Eigen::VectorXd statesConfig = Eigen::VectorXd::Zero(STATES_SIZE,1); //M in the paper

    nh.getParam("states_config", statesConfigXml);
    ROS_ASSERT(statesConfigXml.getType() == XmlRpc::XmlRpcValue::TypeArray);

    if(statesConfigXml.size() != STATES_SIZE)
    {
      ROS_WARN_STREAM("Configuration vector for states_config should have 6 entries.");
    }

    for(int i = 0; i < statesConfigXml.size(); i++)
    {
      statesConfig(i) = static_cast<int>(static_cast<bool>(statesConfigXml[i]));
    }
    return statesConfig;
  }

  double subtractAngles(double angle_k, double angle_k_1)
  {
    angle_k = angle_k*(180/M_PI);
    angle_k_1 = angle_k_1*(180/M_PI);
    double output;

    if((angle_k_1 < 90 && angle_k_1 > 0) && (angle_k > -90 && angle_k < 0))
    {
      output = abs(angle_k) + abs(angle_k_1);
      output = -output;
    }
    else if((angle_k < 90 && angle_k > 0) && (angle_k_1 > -90 && angle_k_1 < 0))
    {
      output = abs(angle_k) + abs(angle_k_1);
    }
    else if((angle_k > 90 && angle_k < 180) && (angle_k_1 < -90 && angle_k_1 > -180))
    {
      output = 360 - (abs(angle_k) + abs(angle_k_1));
      output = -output;
    }
    else if((angle_k_1 > 90 && angle_k_1 < 180) && (angle_k < -90 && angle_k > -180))
    {
      output = 360 - (abs(angle_k) + abs(angle_k_1));
    }
    else
    {
      output = angle_k - angle_k_1;
    }

    return output*(M_PI/180);
  }

  void CovarianceEstimation(int numberOfStates)
  {
    //std::cout << "hamada3" << std::endl << std::endl;
    for(int Idx = 0; Idx < numberOfStates; Idx++)
    {
      for(int Idy = 0; Idy < numberOfStates; Idy++)
      {
        if(ExteroceptiveSensorCovariance_(Idx, Idy) == 0)
          ExteroceptiveSensorCovariance_(Idx, Idy) = 1e-9;
      }
    }

    double MahalabonisConstraint = ProprioceptiveSensorMeasurement_.transpose()
        * ExteroceptiveSensorCovariance_.inverse() *
        ProprioceptiveSensorMeasurement_;

    if(MahalabonisConstraint < MahalanobisThreshold_ * MahalanobisThreshold_)
    {
      ros::Time timeNow = ros::Time::now();
      double timeStamp = timeNow.sec + (timeNow.nsec/1000000000.0);

      static double initialTimeStamp;

      if(firstIteration)
      {
        initialTimeStamp = timeNow.sec + (timeNow.nsec/1000000000.0);
        firstIteration = false;
      }

      timeStamp = timeStamp - initialTimeStamp;
      //std::cout << timeStamp << std::endl << std::endl;

      Eigen::MatrixXd distanceDependentScaleFactorSigmaPoints = Eigen::MatrixXd::Zero(2,(1+2*numberOfStates)*numberOfStates);
      Eigen::MatrixXd timeStampBuffer = Eigen::MatrixXd::Zero(sizeOfTimeBuffer, 2);
      timeStampBuffer.col(0) = Eigen::MatrixXd::Ones(sizeOfTimeBuffer, 1);

      static Eigen::VectorXd previousProprioMeasurement = Eigen::VectorXd::Zero(numberOfStates);
  //    static Eigen::VectorXd previousExteroMeasurement = Eigen::VectorXd::Zero(numberOfStates);

      static Eigen::MatrixXd InnovationMemoryMatrix = Eigen::MatrixXd::Zero(sizeOfTimeBuffer, (2*numberOfStates+1)*numberOfStates);

      static int measurementNumber = 0;

      Eigen::MatrixXd ExteroceptiveSensorCovarianceSqrt = Eigen::MatrixXd::Zero(numberOfStates, numberOfStates);
      Eigen::MatrixXd ExteroSigmaPoints = Eigen::MatrixXd::Zero(numberOfStates, (2*numberOfStates)+1);

      ExteroceptiveSensorCovarianceSqrt = ExteroCovConfidenceFactor_*ExteroceptiveSensorCovariance_.sqrt();
      ExteroSigmaPoints.col(0) = ExteroceptiveSensorMeasurement_;

      for(int i = 1; i < (1+numberOfStates); i++)
        ExteroSigmaPoints.col(i) = ExteroSigmaPoints.col(0) + ExteroceptiveSensorCovarianceSqrt.col(i-1);

      for(int i = 1+numberOfStates; i < 1+(2*numberOfStates); i++)
        ExteroSigmaPoints.col(i) = ExteroSigmaPoints.col(0) - ExteroceptiveSensorCovarianceSqrt.col((i-numberOfStates)-1);

      Eigen::MatrixXd innovationSigmaPoints = Eigen::MatrixXd::Zero(numberOfStates, 2*numberOfStates+1);

      for(int i = 0; i< 2*numberOfStates + 1; i++)
      {
        innovationSigmaPoints.col(i) = ProprioceptiveSensorMeasurement_ - ExteroSigmaPoints.col(i);

        if((StatesConfig_(5) == 1 && StatesConfig_(4) == 0 && StatesConfig_(3) == 0)
           || (StatesConfig_(5) == 0 && StatesConfig_(4) == 1 && StatesConfig_(3) == 0)
           || (StatesConfig_(5) == 0 && StatesConfig_(4) == 0 && StatesConfig_(3) == 1))
        {
          innovationSigmaPoints(numberOfStates-1, i) = subtractAngles(ProprioceptiveSensorMeasurement_(numberOfStates-1),
                                                                      ExteroSigmaPoints(numberOfStates-1,i));

        }
        else if((StatesConfig_(5) == 1 && StatesConfig_(4) == 1 && StatesConfig_(3) == 0)
                || (StatesConfig_(5) == 0 && StatesConfig_(4) == 1 && StatesConfig_(3) == 1)
                || (StatesConfig_(5) == 1 && StatesConfig_(4) == 0 && StatesConfig_(3) == 1))
        {
          innovationSigmaPoints(numberOfStates-1, i) = subtractAngles(ProprioceptiveSensorMeasurement_(numberOfStates-1),
                                                                      ExteroSigmaPoints(numberOfStates-1,i));
          innovationSigmaPoints(numberOfStates-2, i) = subtractAngles(ProprioceptiveSensorMeasurement_(numberOfStates-2),
                                                                      ExteroSigmaPoints(numberOfStates-2,i));
        }
        else if((StatesConfig_(5) == 1 && StatesConfig_(4) == 1 && StatesConfig_(3) == 0))
        {
          innovationSigmaPoints(numberOfStates-1, i) = subtractAngles(ProprioceptiveSensorMeasurement_(numberOfStates-1),
                                                                      ExteroSigmaPoints(numberOfStates-1,i));
          innovationSigmaPoints(numberOfStates-2, i) = subtractAngles(ProprioceptiveSensorMeasurement_(numberOfStates-2),
                                                                      ExteroSigmaPoints(numberOfStates-2,i));
          innovationSigmaPoints(numberOfStates-3, i) = subtractAngles(ProprioceptiveSensorMeasurement_(numberOfStates-3),
                                                                      ExteroSigmaPoints(numberOfStates-3,i));
        }
      }
      //std::cout << innovationSigmaPoints.col(0) << std::endl << std::endl;
      for(int i = 0; i < numberOfStates; i++)
      {
        InnovationMemoryMatrix.block(measurementNumber, i*((2*numberOfStates)+1), 1, (2*numberOfStates + 1))
            = innovationSigmaPoints.row(i);
      }
      timeStampBuffer(measurementNumber,1) = timeStamp;
      measurementNumber++;
      if(measurementNumber == sizeOfTimeBuffer)
      {
        measurementNumber = 0;
      }
      //std::cout << "hamada4" << std::endl << std::endl;
      Eigen::MatrixXd IKtIK = timeStampBuffer.transpose()*timeStampBuffer;

      for(int i = 0; i < (1+2*numberOfStates)*numberOfStates; i++)
      {
        distanceDependentScaleFactorSigmaPoints.col(i) = (IKtIK.inverse()*timeStampBuffer.transpose())*InnovationMemoryMatrix.col(i);
      }

      Eigen::VectorXd distanceDependentScaleFactorMeanWeights = Eigen::VectorXd::Zero((1+2*numberOfStates)*numberOfStates);
      for(int i = 0; i < (1+2*numberOfStates)*numberOfStates; i = i + (1+2*numberOfStates))
      {
        distanceDependentScaleFactorMeanWeights(i) = float(numberOfStates)/(2*numberOfStates + 1);
      }

      for(int i = 0; i < (1+2*numberOfStates)*numberOfStates; i++)
      {
        if(distanceDependentScaleFactorMeanWeights(i) == 0)
        {
          distanceDependentScaleFactorMeanWeights(i) = float(numberOfStates+1)/(4*numberOfStates*numberOfStates + 2*numberOfStates);
        }
      }

      Eigen::VectorXd distanceDependentScaleFactors = Eigen::VectorXd::Zero(numberOfStates);

      for(int i = 0; i < numberOfStates; i++)
      {
        for(int j = i*(1+2*numberOfStates); j < (i+1)*(1+2*numberOfStates); j++)
        {
          distanceDependentScaleFactors(i) = distanceDependentScaleFactors(i) +
              distanceDependentScaleFactorMeanWeights(j) *
              distanceDependentScaleFactorSigmaPoints(1,j);
        }
      }

      for(int i = 0; i < numberOfStates; i++)
      {
        if(isnan(distanceDependentScaleFactors(i)))
          distanceDependentScaleFactors(i) = 1e-9;
      }
      //std::cout << distanceDependentScaleFactors << std::endl << std::endl;

      Eigen::VectorXd deltaProprioCeptiveReading = ProprioceptiveSensorMeasurement_ - previousProprioMeasurement;
      //std::cout << deltaProprioCeptiveReading << std::endl << std::endl;
      previousProprioMeasurement = ProprioceptiveSensorMeasurement_;
      distanceDependentScaleFactors = distanceDependentScaleFactors.cwiseProduct(deltaProprioCeptiveReading);

      scaleFactors_ = distanceDependentScaleFactors;
      //std::cout << distanceDependentScaleFactors << std::endl << std::endl;
    }
  }

  void publishing(double numberOfStates)
  {
    //std::cout << "hamada5" << std::endl << std::endl;
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(numberOfStates,numberOfStates);
    for(int i = 0; i < numberOfStates; i++)
    {
      Q(i,i) = pow(C*scaleFactors_(i),2);
    }
    //std::cout << Q << std::endl << std::endl;

    int checkCounter = 0;
    for(int i = 0; i < numberOfStates; i++)
    {
      if(ProprioceptiveSensorCovariance_(i,i) == 0)
      checkCounter++;
    }
    static Eigen::MatrixXd previous_P = Eigen::MatrixXd::Zero(numberOfStates, numberOfStates);

    if(checkCounter == numberOfStates)
    {
      ProprioceptiveSensorCovariance_ = Q + previous_P;
      //std::cout << previous_P << std::endl << std::endl;
    }
    else
    {
      ProprioceptiveSensorCovariance_ = ProprioceptiveSensorCovariance_ + Q;
    }
    previous_P = ProprioceptiveSensorCovariance_;

    std::cout << Q << std::endl << std::endl;
    nav_msgs::Odometry ProprioCeptiveSensorWithCovariance;
    int Idx;
    int Idy;
    for(int i = 0, Idx = 0; i < STATES_SIZE; i++)
    {
      if(StatesConfig_(i) == 1)
      {
        for(int j = 0, Idy = 0; j < STATES_SIZE; j++)
        {
          if(StatesConfig_(j) == 1)
          {
            ProprioCeptiveSensorWithCovariance.pose.covariance[(i*6)+j] = ProprioceptiveSensorCovariance_(Idx,Idy);
            Idy++;
          }
        }
        Idx++;
      }
    }
    Eigen::VectorXd States_ = Eigen::VectorXd::Zero(6,1);
    for(int i = 0, Idx = 0; i < STATES_SIZE; i++)
    {
      if(StatesConfig_(i) == 1)
      {
        States_(i) = ProprioceptiveSensorMeasurement_(Idx);
        Idx++;
      }
    }
    ProprioCeptiveSensorWithCovariance.pose.pose.position.x = States_(0);
    ProprioCeptiveSensorWithCovariance.pose.pose.position.y = States_(1);
    ProprioCeptiveSensorWithCovariance.pose.pose.position.z = States_(2);

    tf::Quaternion q;
    q.setEuler(States_(3), States_(4), States_(5));

    ProprioCeptiveSensorWithCovariance.pose.pose.orientation.x = q.getX();
    ProprioCeptiveSensorWithCovariance.pose.pose.orientation.y = q.getY();
    ProprioCeptiveSensorWithCovariance.pose.pose.orientation.z = q.getZ();
    ProprioCeptiveSensorWithCovariance.pose.pose.orientation.w = q.getW();

    ProprioCeptiveSensorWithCovariance.header.frame_id = vehicle + "/" + local_frame_id;
    ProprioCeptiveSensorWithCovariance.header.stamp = readingTime;

    odomWithCalculatedCovPublisher.publish(ProprioCeptiveSensorWithCovariance);
    //std::cout << "hamada6" << std::endl << std::endl;
  }

  void odomWithNoCovCb(const nav_msgs::OdometryConstPtr &odomWithNoCovMsg)
  {
    //std::cout << "hamada1" << std::endl << std::endl;
    int NumberOfStates = StatesConfig_.sum();
    ProprioceptiveSensorMeasurement_ = Eigen::VectorXd::Zero(NumberOfStates);
    ProprioceptiveSensorCovariance_ = Eigen::MatrixXd::Zero(NumberOfStates,NumberOfStates);

    Eigen::Matrix <double, STATES_SIZE, 1> odometryMessagePose;
    readingTime = odomWithNoCovMsg->header.stamp;

    odometryMessagePose(0) = odomWithNoCovMsg->pose.pose.position.x;
    odometryMessagePose(1) = odomWithNoCovMsg->pose.pose.position.y;
    odometryMessagePose(2) = odomWithNoCovMsg->pose.pose.position.z;

    Eigen::Matrix <double, STATES_SIZE/2, 1> orientation;
    tf::Quaternion q(odomWithNoCovMsg->pose.pose.orientation.x, odomWithNoCovMsg->pose.pose.orientation.y,
                     odomWithNoCovMsg->pose.pose.orientation.z, odomWithNoCovMsg->pose.pose.orientation.w);
    tf::Matrix3x3 m(q);
    m.getRPY(odometryMessagePose(3),odometryMessagePose(4),odometryMessagePose(5));

    int Idx = 0;
    int IdxMeasurement = 0;
    while(Idx < STATES_SIZE)
    {
      if(StatesConfig_(Idx) == 1)
      {
        ProprioceptiveSensorMeasurement_(IdxMeasurement) = odometryMessagePose(Idx);
        IdxMeasurement++;
        Idx++;
      }
      else
      {
        Idx++;
      }
    }

    int IdxxMeasurement = 0;
    int IdxyMeasurement = 0;
    for(int Idxx = 0, IdxxMeasurement = 0; Idxx < STATES_SIZE; Idxx++)
    {
      if(StatesConfig_(Idxx) == 1)
      {
        for(int Idxy = 0, IdxyMeasurement = 0; Idxy < STATES_SIZE; Idxy++)
        {
          if(StatesConfig_(Idxy) == 1)
          {
            ProprioceptiveSensorCovariance_(IdxxMeasurement, IdxyMeasurement) = odomWithNoCovMsg->pose.covariance[(Idxx*STATES_SIZE)+Idxy];
            IdxyMeasurement++;
          }
        }
        IdxxMeasurement++;
      }
    }
    //std::cout << "hamada2" << std::endl << std::endl;
    publishing(NumberOfStates);
  }

  void odomWithCovCb(const nav_msgs::OdometryConstPtr &odomWithCovMsg)
  {
    int NumberOfStates = StatesConfig_.sum();
    ExteroceptiveSensorMeasurement_ = Eigen::VectorXd::Zero(NumberOfStates);
    ExteroceptiveSensorCovariance_ = Eigen::MatrixXd::Zero(NumberOfStates,NumberOfStates);

    Eigen::Matrix <double, STATES_SIZE, 1> odometryMessagePose;
    odometryMessagePose(0) = odomWithCovMsg->pose.pose.position.x;
    odometryMessagePose(1) = odomWithCovMsg->pose.pose.position.y;
    odometryMessagePose(2) = odomWithCovMsg->pose.pose.position.z;

    Eigen::Matrix <double, STATES_SIZE/2, 1> orientation;
    tf::Quaternion q(odomWithCovMsg->pose.pose.orientation.x, odomWithCovMsg->pose.pose.orientation.y,
                     odomWithCovMsg->pose.pose.orientation.z, odomWithCovMsg->pose.pose.orientation.w);
    tf::Matrix3x3 m(q);
    m.getRPY(odometryMessagePose(3),odometryMessagePose(4),odometryMessagePose(5));

    int Idx = 0;
    int IdxMeasurement = 0;
    while(Idx < STATES_SIZE)
    {
      if(StatesConfig_(Idx) == 1)
      {
        ExteroceptiveSensorMeasurement_(IdxMeasurement) = odometryMessagePose(Idx);
        IdxMeasurement++;
        Idx++;
      }
      else
      {
        Idx++;
      }
    }

    int IdxxMeasurement = 0;
    int IdxyMeasurement = 0;
    for(int Idxx = 0, IdxxMeasurement = 0; Idxx < STATES_SIZE; Idxx++)
    {
      if(StatesConfig_(Idxx) == 1)
      {
        for(int Idxy = 0, IdxyMeasurement = 0; Idxy < STATES_SIZE; Idxy++)
        {
          if(StatesConfig_(Idxy) == 1)
          {
            ExteroceptiveSensorCovariance_(IdxxMeasurement, IdxyMeasurement) = odomWithCovMsg->pose.covariance[(Idxx*STATES_SIZE)+Idxy];
            IdxyMeasurement++;
          }
        }
        IdxxMeasurement++;
      }
    }
    CovarianceEstimation(NumberOfStates);
  }

  int main(int argc, char **argv)
  {
    ros::init(argc, argv, "covariance_estimation");
    ros::NodeHandle nh, nh_private("~");

    std::string odomWithNoCovTopic;
    std::string odomWithCovTopic;

    nh_private.getParam("odom_with_no_covariance_topic", odomWithNoCovTopic);
    nh_private.getParam("odom_with_covariance_topic", odomWithCovTopic);

    ros::Subscriber odomWithNoCovSub = nh.subscribe(odomWithNoCovTopic, 10, odomWithNoCovCb);
    ros::Subscriber odomWithCovSub = nh.subscribe(odomWithCovTopic, 1, odomWithCovCb);

    std::string odomWithCalculatedCovTopic;

    nh_private.getParam("odom_with_calculated_cov_topic", odomWithCalculatedCovTopic);

    odomWithCalculatedCovPublisher = nh.advertise<nav_msgs::Odometry>(odomWithCalculatedCovTopic, 10);

    StatesConfig_ = loadStatesConfiguration(nh_private);

    //reading the algorithm parameters
    nh_private.getParam("size_of_Xi", sizeOfTimeBuffer);
    nh_private.getParam("cov_confidence_factor", ExteroCovConfidenceFactor_);
    nh_private.getParam("Mahalanobis_threshold", MahalanobisThreshold_);
    nh_private.getParam("standard_deviation_correlation_constant_distance", C);
    nh.getParam("vehicle", vehicle);
    nh.getParam("local_frame_id", local_frame_id);

    //std::cout << StatesConfig_ << std::endl << std::endl;
    std::cout << StatesConfig_ << std::endl << std::endl;

    int numberOfStates = StatesConfig_.sum();

    scaleFactors_ = Eigen::VectorXd::Zero(numberOfStates);
    ExteroceptiveSensorMeasurement_ = Eigen::VectorXd::Zero(numberOfStates);
    ExteroceptiveSensorCovariance_ = Eigen::MatrixXd::Zero(numberOfStates, numberOfStates);

    ros::spin();

    return EXIT_SUCCESS;
  }
