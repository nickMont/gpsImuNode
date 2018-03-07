//BORKED DO NOT USE

#include "filter.h"
#include <Eigen/LU>  // For matrix inverse
#include "fdesSetter.h"
#include "ros/ros.h"
#include "sensor_msgs/Image.h"

namespace gps_odom {


class storeControlInfo
{
public:
  storeControlInfo(std::string filename)
  {
    std::string quadPoseTopic, quadName;
    quadName = ros::this_node::getName();
    sub_ = n_.subscribe(quadName+"/Fdes", 1, &storeControlInfo::callback, this);
  }
  void callback(const sensor_msgs::Image& imgmsg)
  {
    Fdes_t=Fdes_t+1;
  }

  private:
  ros::NodeHandle n_; 
  ros::Subscriber sub_;

};//End of class

int main(int argc, char **argv)
{

  ros::init(argc, argv, "storeControlInfo");
  storeControlInfo controlStorage("derp");
  ros::spin();

  return 0;
}


}  // namespace gps_odom
