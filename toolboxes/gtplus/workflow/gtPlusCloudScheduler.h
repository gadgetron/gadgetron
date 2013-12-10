/** \file   gtPlusCloudScheduler.h
    \brief  Define and implement the GtPlus cloud job scheduler class
            A simple scheduling strategy is implemented here. The number of job packages which are sent
            to a node is propotional to the computing power index for that node.

            This class may serve as the base class to implement more complicated job scheduling strategies.

    \author Hui Xue
*/

#pragma once

#include "GtPlusExport.h"
#include "gtPlusISMRMRDReconUtil.h"

namespace Gadgetron { namespace gtPlus {

/**
The scheduler class for gadgetron cloud.
This class can serves as the base class for more complicated scheduling strategy.
*/

class EXPORTGTPLUS gtPlusCloudScheduler
{
public:

    gtPlusCloudScheduler();
    virtual ~gtPlusCloudScheduler();

    virtual void printInfo(std::ostream& os) const;

    // compute the scheduling for every job
    // nodeIDforJobs stores the node ID to run every job
    // node ID starts from 0
    virtual bool schedulerJobs(std::vector<int>& nodeIDforJobs);

    void setNumOfJobs(size_t numOfJobs);

    void setUpNodes(size_t numOfNodes);
    void setUpNodes(const std::vector<double>& nodeComputingPowerIndexes);

protected:

    // number of nodes
    size_t num_of_nodes_;

    // number of jobs, for this simple scheduler, all jobs are considered to have equal sizes
    size_t num_of_jobs_;

    // computing power indexes for every nodes; if not set, all nodes are treated to have equal computing powers
    std::vector<std::pair<int, double> > node_id_computing_power_indexes_;
};

}}
