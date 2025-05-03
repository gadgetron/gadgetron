#include "StandardDeviationGadget.h"

using namespace Gadgetron::Indexing;
namespace Gadgetron {

StandardDeviationGadget::StandardDeviationGadget(const Core::Context& context, const Core::GadgetProperties& props)
    : ChannelGadget(context, props) {
    header = context.header;
    image_counter_ = 0;
}

void StandardDeviationGadget::process(Core::InputChannel<IsmrmrdReconData>& input, Core::OutputChannel& out) {
    std::vector<IsmrmrdReconData> pseudoreplicaList;
    for (IsmrmrdReconData reconData : input) {
        pseudoreplicaList.push_back(reconData);
    }

    GDEBUG("ReconDataList Length: %d\n", pseudoreplicaList.size());
    IsmrmrdDataBuffered& dbuff = pseudoreplicaList[0].rbit_.begin()->data_;
    uint16_t REP = pseudoreplicaList.size();
    uint16_t E0 = dbuff.data_.get_size(0);
    uint16_t E1 = dbuff.data_.get_size(1);
    uint16_t E2 = dbuff.data_.get_size(2);
    uint16_t CHA = dbuff.data_.get_size(3);
    uint16_t N = dbuff.data_.get_size(4);
    uint16_t S = dbuff.data_.get_size(5);
    uint16_t LOC = dbuff.data_.get_size(6);
    GDEBUG("%d,%d,%d,%d,%d,%d,%d,%d \n", REP, E0, E1, E2, CHA, N, S, LOC);

    hoNDArray<std::complex<float>> standardDeviations(E0, E1, E2, CHA, N, S, LOC);
    for (size_t e0 = 0; e0 < E0; e0++) {
        for (size_t e1 = 0; e1 < E1; e1++) {
            for (size_t e2 = 0; e2 < E2; e2++) {
                for (size_t cha = 0; cha < CHA; cha++) {
                    for (size_t n = 0; n < N; n++) {
                        for (size_t s = 0; s < S; s++) {
                            for (size_t loc = 0; loc < LOC; loc++) {
                                hoNDArray<std::complex<float>> replicaData(REP);
                                for (size_t rep = 0; rep < REP; rep++) {
                                    replicaData[rep] =
                                        pseudoreplicaList[rep].rbit_[0].data_.data_[e0, e1, e2, cha, n, s, loc];
                                }
                                standardDeviations[e0, e1, e2, cha, n, s, loc] = stddev(replicaData);
                            }
                        }
                    }
                }
            }
        }
    }

    auto meanOfStdDevs = mean(abs(standardDeviations));
    GDEBUG("Mean of Std. Deviations: %f\n", meanOfStdDevs);
    if(meanOfStdDevs > errorThreshold){
        GERROR("Mean of Standard Deviations exceeds threshold: %f > %f",meanOfStdDevs, errorThreshold);
    }
    pseudoreplicaList[0].rbit_[0].data_.data_ = standardDeviations;
    out.push(std::move(pseudoreplicaList[0]));
}

GADGETRON_GADGET_EXPORT(StandardDeviationGadget);

} // namespace Gadgetron
