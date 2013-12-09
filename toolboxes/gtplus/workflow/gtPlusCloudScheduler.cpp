/** \file   gtPlusCloudScheduler.cpp
    \brief  Define and implement the GtPlus cloud job scheduler class
            A simple scheduling strategy is implemented here. The number of job packages which are sent
            to a node is propotional to the computing power index for that node.

            This class may serve as the base class to implement more complicated job scheduling strategies.

    \author Hui Xue
*/

#include "gtPlusCloudScheduler.h"

namespace Gadgetron { namespace gtPlus {

gtPlusCloudScheduler::gtPlusCloudScheduler() : num_of_nodes_(0), num_of_jobs_(0)
{
}

gtPlusCloudScheduler::~gtPlusCloudScheduler()
{
}

void gtPlusCloudScheduler::printInfo(std::ostream& os) const
{
    using namespace std;

    os << "-------------- GTPlus Cloud scheduler for jobs ---------------" << endl;
    os << "This class implements the simple scheduling scheme for GtPlus cloud " << endl;
    os << "The scheduler here tries to allocate nodes to jobs propotional to the power indexes provided " << endl;
    os << "--------------------------------------------------------------" << endl;
}

void gtPlusCloudScheduler::setNumOfJobs(unsigned long long numOfJobs)
{
    num_of_jobs_ = numOfJobs;
}

void gtPlusCloudScheduler::setUpNodes(unsigned long long numOfNodes)
{
    num_of_nodes_ = numOfNodes;
    if ( num_of_nodes_ > 0 )
    {
        node_id_computing_power_indexes_.resize(num_of_nodes_);
        for ( unsigned long long ii=0; ii<num_of_nodes_; ii++ )
        {
            node_id_computing_power_indexes_[ii].first = ii;
            node_id_computing_power_indexes_[ii].second = 1.0;
        }
    }
}

void gtPlusCloudScheduler::setUpNodes(const std::vector<double>& nodeComputingPowerIndexes)
{
    num_of_nodes_ = nodeComputingPowerIndexes.size();
    node_id_computing_power_indexes_.resize(num_of_nodes_);

    for ( unsigned long long ii=0; ii<num_of_nodes_; ii++ )
    {
        node_id_computing_power_indexes_[ii].first = (int)ii;
        node_id_computing_power_indexes_[ii].second = nodeComputingPowerIndexes[ii];
    }
}

struct gtPlusCloudSchedulerNodeSorter
{
    gtPlusCloudSchedulerNodeSorter() {}
    ~gtPlusCloudSchedulerNodeSorter() {}

    bool operator()(const std::pair<int, double>& A, const std::pair<int, double>& B) const
    {
        return (A.second > B.second);
    }
};

bool gtPlusCloudScheduler::schedulerJobs(std::vector<int>& nodeIDforJobs)
{
    try
    {
        unsigned long long ii;

        nodeIDforJobs.clear();

        if ( num_of_nodes_==0 || num_of_jobs_==0 )
        {
            GADGET_WARN_MSG("num_of_nodes_==0 || num_of_jobs_==0");
            return true;
        }

        if ( node_id_computing_power_indexes_.size() < num_of_nodes_ )
        {
            GADGET_WARN_MSG("node_computing_power_indexes_.size() < num_of_nodes_ : computing power indexes for all nodes are set to be equal ... ");
            node_id_computing_power_indexes_.resize(num_of_nodes_, std::pair<int, double>(0, 1.0) );
            for ( ii=0; ii<num_of_nodes_; ii++ )
            {
                node_id_computing_power_indexes_[ii].first = (int)ii;
            }
        }

        nodeIDforJobs.resize(num_of_jobs_, -1);

        // always sort the nodes with higher computing power node ahead
        std::sort(node_id_computing_power_indexes_.begin(), node_id_computing_power_indexes_.end(), gtPlusCloudSchedulerNodeSorter() );

        if ( num_of_jobs_ <= num_of_nodes_ )
        {
            for ( ii=0; ii<num_of_jobs_; ii++ )
            {
                nodeIDforJobs[ii] = node_id_computing_power_indexes_[ii].first;
            }
        }
        else
        {
            double totalComputingPower = 0.0;
            for ( ii=0; ii<num_of_nodes_; ii++ )
            {
                totalComputingPower += node_id_computing_power_indexes_[ii].second;
            }

            unsigned long long totalJobAllocated = 0;
            std::vector<unsigned long long> jobPerNode(num_of_nodes_, 0);
            for ( ii=0; ii<num_of_nodes_; ii++ )
            {
                jobPerNode[ii] = (unsigned long long)(std::floor(num_of_jobs_ * node_id_computing_power_indexes_[ii].second/totalComputingPower));
                totalJobAllocated += jobPerNode[ii];
            }

            if ( totalJobAllocated < num_of_jobs_ )
            {
                // give high computing power nodes more jobs
                for ( ii=0; ii<(num_of_jobs_-totalJobAllocated); ii++ )
                {
                    jobPerNode[ii%num_of_nodes_]++;
                }
            }

            unsigned long long jobID = 0;
            for ( ii=0; ii<num_of_nodes_; ii++ )
            {
                for ( unsigned long long jj=0; jj<jobPerNode[ii]; jj++ )
                {
                    nodeIDforJobs[jobID++] = node_id_computing_power_indexes_[ii].first;
                }
            }

            GADGET_CHECK_RETURN_FALSE(jobID==num_of_jobs_);
        }
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusCloudScheduler::schedulerJobs(std::vector<int>& nodeIDforJobs) ... ");
        return false;
    }

    return true;
}

}}
