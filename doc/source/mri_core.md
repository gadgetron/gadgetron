# MRI Core Gadgets
The MRI core gadgets are the backbone of the Gadgetron pipeline, and provide core functions like data-ingress, flow management, type conversion, and more. 

## Pure Gadgets
| Gadget Name                         | Has Tests? | Inputs                                           |
| ----------------------------------- | ---------- | ------------------------------------------------ |
| AutoScaleGadget                     | Yes        | ImageHeader and hoNDArray (float)                |
| ComplexToFloatGadget                | Yes        | Image (complex float)                            |
| DenoiseGadget                       | Yes        | DenoiseImage (float/complex float or ImageArray) |
| ScaleGadget                         | No         | Image (float) or IsmrmrdImageArray               |

## Channel Gadgets
| Gadget Name                         | Has Tests? | Inputs                                           |
| ----------------------------------- | ---------- | ------------------------------------------------ |
| AcquisitionAccumulateTriggerGadget  | Yes        | Acquisition or Waveform                          |
| AcquisitionPassthroughGadget        | Yes        | Acquisition                                      |
| AsymmetricEchoAdjustROGadget        | Yes        | Acquisition                                      |
| BucketToBufferGadget                | Yes        | AcquisitionBucket                                |
| CplxDumpGadget                      | No         | Acquisition                                      |
| ExtractGadget                       | Yes        | Image (complex float)                            |
| FlagTriggerGadget                   | No         | Acquisition                                      |
| FloatToFixPointGadget               | No         | Image (float)                                    |
| ImageFinishGadget                   | Yes        | ------------------------------------------------ |
| ImageIndexGadget                    | Yes        | Image (any)                                      |
| ImageResizingGadget                 | No         | Image (any)                                      |
| IsmrmrdDumpGadget                   | No         | Acquisition or Waveform                          |
| NoiseAdjustGadget                   | Yes        | Acquisition                                      |
| PhysioInterpolationGadget           | No         | Image (complex float)                            |

## Gadget1 Gadgets
| Gadget Name                         | Has Tests? | Inputs                                           |
| ----------------------------------- | ---------- | ------------------------------------------------ |
| FFTGadget                           | No         | IsmrmrdReconData                                 |
| ImageAccumulatorGadget              | Yes        | IsmrmrdImageArray                                |
| ImageSortGadget                     | Yes        | ImageHeader                                      |
| PseudoReplicatorGadget              | No         | IsmrmrdReconData                                 |
| SimpleReconGadget                   | Yes        | IsmrmrdReconData                                 |

## Gadget2 Gadgets
| Gadget Name                         | Has Tests? | Inputs                                           |
| ----------------------------------- | ---------- | ------------------------------------------------ |
| AccumulatorGadget                   | No         | AcquisitionHeader and hoNDArray (complex float)  |
| CoilComputationGadget               | No         | ImageHeader and hoNDArray (complex float)        |
| CoilReductionGadget                 | Yes        | AcquisitionHeader and hoNDArray (complex float)  |
| CombineGadget                       | No         | ImageHeader and hoNDArray (complex float)        |
| CropAndCombineGadget                | No         | ImageHeader and hoNDArray (complex float)        |
| FlowPhaseSubtractionGadget          | No         | ImageHeader and hoNDArray (complex float)        |
| ImageWriterGadget                   | Yes        | ImageHeader and hoNDArray (short/float/complex)  |
| MaxwellCorrectionGadget             | No         | ImageHeader and hoNDArray (complex float)        |
| NoiseAdjustGadget_unoptimized       | No         | AcquisitionHeader and hoNDArray (complex float)  |
| PCACoilGadget                       | Yes        | AcquisitionHeader and hoNDArray (complex float)  |
| PartialFourierAdjustROGadget        | No         | AcquisitionHeader and hoNDArray (complex float)  |
| RemoveROOversamplingGadget          | Yes        | AcquisitionHeader and hoNDArray (complex float)  |
| WhiteNoiseInjectorGadget            | No         | AcquisitionHeader and hoNDArray (complex float)  |

## Gadget3 Gadgets
| Gadget Name                         | Has Tests? | Inputs                                           |
| ----------------------------------- | ---------- | ------------------------------------------------ |
| ImageFFTGadget                      | No         | ImageHeader, hoNDArray, and MetaContainer        |

## Other Gadgets
| Gadget Name            | Type                | Has Tests? | Inputs                                           |
| ---------------------- | ------------------- | ---------- | ------------------------------------------------ |
| ImageArraySplitGadget  | Gadget1Of2          | No         | ImageHeader and IsmrmrdImageArray                |
| RateLimitGadget        | BasicPropertyGadget | No         | ACEMessageBlock                                  |

## Generic Recon Gadgets
| Generic Recon Gadgets                          | Has Tests? | Inputs                                           |
| ---------------------------------------------- | ---------- | ------------------------------------------------ |
| GenericImageReconArrayToImageGadget            | No         | ------------------------------------------------ |
| GenericImageReconGadget                        | No         | ------------------------------------------------ |
| GenericReconAccumulateImageTriggerGadget       | No         | ------------------------------------------------ |
| GenericReconCartesianFFTGadget                 | No         | ------------------------------------------------ |
| GenericReconCartesianGrappaAIGadget            | No         | ------------------------------------------------ |
| GenericReconCartesianGrappaGadget              | Yes        | ------------------------------------------------ |
| GenericReconCartesianNonLinearSpirit2DTGadget  | No         | ------------------------------------------------ |
| GenericReconCartesianReferencePrepGadget       | Yes        | ------------------------------------------------ |
| GenericReconCartesianSpiritGadget              | No         | ------------------------------------------------ |
| GenericReconEigenChannelGadget                 | Yes        | ------------------------------------------------ |
| GenericReconFieldOfViewAdjustmentGadget        | Yes        | ------------------------------------------------ |
| GenericReconGadget                             | No         | ------------------------------------------------ |
| GenericReconImageArrayScalingGadget            | Yes        | ------------------------------------------------ |
| GenericReconImageToImageArrayGadget            | No         | ------------------------------------------------ |
| GenericReconKSpaceFilteringGadget              | Yes        | ------------------------------------------------ |
| GenericReconNoiseStdMapComputingGadget         | No         | ------------------------------------------------ |
| GenericReconPartialFourierHandlingFilterGadget | No         | ------------------------------------------------ |
| GenericReconPartialFourierHandlingGadget       | No         | ------------------------------------------------ |
| GenericReconPartialFourierHandlingPOCSGadget   | Yes        | ------------------------------------------------ |
| GenericReconReferenceKSpaceDelayedBufferGadget | No         | ------------------------------------------------ |
