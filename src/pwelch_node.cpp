#include <ros/ros.h>
#include "rosneuro_pwelch/PwelchNode.h"
#include "rosneuro_buffers_ringbuffer/RingBuffer.hpp"

int main(int argc, char** argv) {

	ros::init(argc, argv, "pwelch_node");

	rosneuro::PwelchNode<double> pwelchNode;

	if(pwelchNode.configure() == false) {
		ROS_ERROR("[Pwelch] Configuration failed");
		return -1;
	}

	pwelchNode.run();

	ros::shutdown();


	return 0;
}
