#include "ace/INET_Addr.h"
#include "ace/SOCK_Stream.h"
#include "ace/SOCK_Connector.h"
#include "ace/Log_Msg.h"
#include "ace/Get_Opt.h"
#include "ace/OS_NS_string.h"
#include "ace/Reactor.h"

#include "../gadgetheaders.h"
#include "siemensraw.hpp"
#include "GadgetSocketReceiver.h"
#include "ImageWriter.h"
#include "ConfigParser.h"
#include "NDArray.h"

struct spiral_parameters {
  int Interleaves;
  int ADCsPerInterleave;
  int SamplesPerADC;
  int SamplesToSkipStart;
  int SamplesToSkipEnd; 
  int SamplingTime_ns;
  int Reordering;
  double MaxGradient_Gcm;
  double MaxSlewRate_Gcms;
  double krmax_cm;
  double FOVCoeff_1;
};

int getSpiralParameters(SiemensRawData* sd, spiral_parameters* parm)
{
  std::string val;
  if (sd->GetMeasYapsParameter(std::string("sKSpace.lRadialViews"), val) == -1) {
    ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Unable to findsKSpace.lRadialViews\n")) );
    return -1;
  }
  parm->Interleaves = atoi(val.c_str());

  if (sd->GetMeasYapsParameter(std::string("sWiPMemBlock.alFree[52]"), val) == -1) {
    ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Unable to find sWiPMemBlock.alFree[52]\n")) );
    return -1;
  }
  parm->ADCsPerInterleave = atoi(val.c_str());

  if (sd->GetMeasYapsParameter(std::string("sWiPMemBlock.alFree[53]"), val) == -1) {
    ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Unable to find sWiPMemBlock.alFree[53]\n")) );
    return -1;
  }
  parm->SamplesPerADC = atoi(val.c_str());

  if (sd->GetMeasYapsParameter(std::string("sWiPMemBlock.alFree[54]"), val) == -1) {
    //ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Unable to find sWiPMemBlock.alFree[54]\n")) );
    parm->SamplesToSkipStart = 0;//atoi(val.c_str());
  } else {
    parm->SamplesToSkipStart = atoi(val.c_str());
  }

  if (sd->GetMeasYapsParameter(std::string("sWiPMemBlock.alFree[55]"), val) == -1) {
    ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Unable to find sWiPMemBlock.alFree[55]\n")) );
    return -1;
  }
  parm->SamplesToSkipEnd = atoi(val.c_str());

  if (sd->GetMeasYapsParameter(std::string("sWiPMemBlock.alFree[56]"), val) == -1) {
    ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Unable to find sWiPMemBlock.alFree[56]\n")) );
    return -1;
  }
  parm->SamplingTime_ns = atoi(val.c_str());

  if (sd->GetMeasYapsParameter(std::string("sKSpace.unReordering"), val) == -1) {
    ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Unable to find sKSpace.unReordering\n")) );
    return -1;
  }
  parm->Reordering = strtol (val.c_str(), NULL, 16);

  if (sd->GetMeasYapsParameter(std::string("sWiPMemBlock.adFree[6]"), val) == -1) {
    ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Unable to find sWiPMemBlock.adFree[6]\n")) );
    return -1;
  }
  parm->MaxGradient_Gcm = atof(val.c_str());

  if (sd->GetMeasYapsParameter(std::string("sWiPMemBlock.adFree[7]"), val) == -1) {
    ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Unable to find sWiPMemBlock.adFree[7]\n")) );
    return -1;
  }
  parm->MaxSlewRate_Gcms = atof(val.c_str());

  if (sd->GetMeasYapsParameter(std::string("sWiPMemBlock.adFree[8]"), val) == -1) {
    ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Unable to find sWiPMemBlock.adFree[8]\n")) );
    return -1;
  }
  parm->krmax_cm= atof(val.c_str());

  if (sd->GetMeasYapsParameter(std::string("sWiPMemBlock.adFree[9]"), val) == -1) {
    ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Unable to find sWiPMemBlock.adFree[9]\n")) );
    return -1;
  }
  parm->FOVCoeff_1= atof(val.c_str());

  return 0;
}

int dumpSpiralParameters(spiral_parameters* spi_parm) {

  ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Interleaves           = %d\n"), spi_parm->Interleaves) );
  ACE_DEBUG( (LM_DEBUG, ACE_TEXT("ADCsPerInterleave     = %d\n"), spi_parm->ADCsPerInterleave) );
  ACE_DEBUG( (LM_DEBUG, ACE_TEXT("SamplesPerADC         = %d\n"), spi_parm->SamplesPerADC) );
  ACE_DEBUG( (LM_DEBUG, ACE_TEXT("SamplesToSkipStart    = %d\n"), spi_parm->SamplesToSkipStart) );
  ACE_DEBUG( (LM_DEBUG, ACE_TEXT("SamplesToSkipEnd      = %d\n"), spi_parm->SamplesToSkipEnd) );
  ACE_DEBUG( (LM_DEBUG, ACE_TEXT("SamplingTime_ns       = %d\n"), spi_parm->SamplingTime_ns) );
  ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Reordering            = %d\n"), spi_parm->Reordering) );
  ACE_DEBUG( (LM_DEBUG, ACE_TEXT("MaxGradient_Gcm       = %f\n"), spi_parm->MaxGradient_Gcm) );
  ACE_DEBUG( (LM_DEBUG, ACE_TEXT("MaxSlewRate_Gcms      = %f\n"), spi_parm->MaxSlewRate_Gcms) );
  ACE_DEBUG( (LM_DEBUG, ACE_TEXT("krmax_cm              = %f\n"), spi_parm->krmax_cm) );
  ACE_DEBUG( (LM_DEBUG, ACE_TEXT("FOVCoeff_1            = %f\n"), spi_parm->FOVCoeff_1) );
  return 0;
}

int ACE_TMAIN(int argc, ACE_TCHAR *argv[] )
{
  static const ACE_TCHAR options[] = ACE_TEXT(":p:h:c:f:l:");
  
  ACE_Get_Opt cmd_opts(argc, argv, options);
  
  ACE_TCHAR port_no[1024];
  ACE_OS_String::strncpy(port_no, "9002", 1024);
  
  ACE_TCHAR hostname[1024];
  ACE_OS_String::strncpy(hostname, "localhost", 1024);
  
  ACE_TCHAR filename[4096];
  ACE_OS_String::strncpy(filename, "./data.dat", 4096);

  ACE_TCHAR config_lib[1024];
  ACE_OS_String::strncpy(config_lib, "core", 1024);

  ACE_TCHAR config_name[1024];
  ACE_OS_String::strncpy(config_name, "default", 1024);
 
   int option;
  while ((option = cmd_opts()) != EOF) {
    switch (option) {
    case 'p':
      ACE_OS_String::strncpy(port_no, cmd_opts.opt_arg(), 1024);
      break;
    case 'h':
      ACE_OS_String::strncpy(hostname, cmd_opts.opt_arg(), 1024);
      break;
    case 'l':
      ACE_OS_String::strncpy(config_lib, cmd_opts.opt_arg(), 1024);
      break;
    case 'c':
      ACE_OS_String::strncpy(config_name, cmd_opts.opt_arg(), 1024);
      break;
    case 'f':
      ACE_OS_String::strncpy(filename, cmd_opts.opt_arg(), 4096);
      break;
    case ':':
      ACE_ERROR_RETURN((LM_ERROR, ACE_TEXT("-%c requires an argument.\n"), cmd_opts.opt_opt()),-1);
      break;
    default:
      ACE_ERROR_RETURN( (LM_ERROR, ACE_TEXT("Command line parse error\n")), -1);
    }    
  }
  
  ACE_DEBUG(( LM_INFO, ACE_TEXT("Gadgetron Data Sender\n") ));
  
  //We need some data to work with
  SiemensRawData sd;
  sd.ReadRawFile(filename);

  SiemensBaseParameters bp = sd.GetBaseParameters();

  GadgetMessageAcquisition acq_head_base;
  memset(&acq_head_base, 0, sizeof(GadgetMessageAcquisition) ); 

  acq_head_base.min_idx.line = sd.GetMinValues()->sLC.ushLine;
  acq_head_base.min_idx.acquisition = sd.GetMinValues()->sLC.ushAcquisition;
  acq_head_base.min_idx.slice = sd.GetMinValues()->sLC.ushSlice;
  acq_head_base.min_idx.partition = sd.GetMinValues()->sLC.ushPartition;
  acq_head_base.min_idx.echo = sd.GetMinValues()->sLC.ushEcho;
  acq_head_base.min_idx.phase = sd.GetMinValues()->sLC.ushPhase;
  acq_head_base.min_idx.repetition = sd.GetMinValues()->sLC.ushRepetition;
  acq_head_base.min_idx.set = sd.GetMinValues()->sLC.ushSet;
  acq_head_base.min_idx.segment = sd.GetMinValues()->sLC.ushSeg;
  acq_head_base.min_idx.channel = sd.GetMinValues()->ushChannelId;

  acq_head_base.max_idx.line = sd.GetMaxValues()->sLC.ushLine;
  acq_head_base.max_idx.acquisition = sd.GetMaxValues()->sLC.ushAcquisition;
  acq_head_base.max_idx.slice = sd.GetMaxValues()->sLC.ushSlice;
  acq_head_base.max_idx.partition = sd.GetMaxValues()->sLC.ushPartition;
  acq_head_base.max_idx.echo = sd.GetMaxValues()->sLC.ushEcho;
  acq_head_base.max_idx.phase = sd.GetMaxValues()->sLC.ushPhase;
  acq_head_base.max_idx.repetition = sd.GetMaxValues()->sLC.ushRepetition;
  acq_head_base.max_idx.set = sd.GetMaxValues()->sLC.ushSet;
  acq_head_base.max_idx.segment = sd.GetMaxValues()->sLC.ushSeg;
  acq_head_base.max_idx.channel = sd.GetMaxValues()->ushChannelId;

  GadgetMessageIdentifier id;
  GadgetMessageConfigurator conf;

  id.id = GADGET_MESSAGE_CONFIGURATION;
  ACE_OS_String::strncpy(conf.configurator_lib,config_lib,1024);
  ACE_OS_String::strncpy(conf.configurator_name,config_name,1024);


  std::string ucTrajectory;
  int trajectory;
  spiral_parameters spi_parm;
  if (sd.GetMeasYapsParameter(std::string("sKSpace.ucTrajectory"), ucTrajectory) == -1) {
    ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Unable to find trajectory parameter, assuming Cartesian data\n")) ); 
    trajectory = TRAJECTORY_CARTESIAN;
  } else {
    trajectory = strtol (ucTrajectory.c_str(), NULL, 16);
    ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Trajectory : %d\n"), trajectory) );
  }


  ConfigParser cp;

  switch (trajectory) {

  case TRAJECTORY_CARTESIAN:
    cp.add(std::string("encoding"), std::string("trajectory"),
	   std::string("cartesian"));
    break;
  case TRAJECTORY_RADIAL:
    cp.add(std::string("encoding"), std::string("trajectory"),
	   std::string("radial"));
    break;
  case TRAJECTORY_SPIRAL:
    cp.add(std::string("encoding"), std::string("trajectory"),
	   std::string("spiral"));
    break;
  default:
    ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Unknown trajectory, assumin Cartesian\n")) );
    cp.add(std::string("encoding"), std::string("trajectory"),
	   std::string("cartesian"));
    break;
  };
  
  cp.add(std::string("encoding"), std::string("matrix_x"), 
	 (size_t)bp.matrix_size[0]);

  cp.add(std::string("encoding"), std::string("matrix_y"), 
	 (size_t)bp.matrix_size[1]);

  cp.add(std::string("encoding"), std::string("matrix_z"), 
	 (size_t)bp.matrix_size[2]);

  cp.add(std::string("encoding"), std::string("readout_length"), 
	 (size_t)sd.GetFirstNode()->mdh.ushSamplesInScan);

  cp.add(std::string("encoding"), std::string("channels"), 
	 (size_t)sd.GetFirstNode()->mdh.ushUsedChannels);

  cp.add(std::string("encoding"), std::string("base_resolution"), 
	 (size_t)bp.base_resolution);

  cp.add(std::string("encoding"), std::string("slices"), 
	 (size_t)acq_head_base.max_idx.slice+1);


  if (trajectory == TRAJECTORY_SPIRAL) {
    if (getSpiralParameters(&sd,&spi_parm) < 0) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("Unable to locate spiral parameters\n")) );
      return -1;
    } else {
      cp.add(std::string("spiral"), std::string("Interleaves"), 
	     (size_t)spi_parm.Interleaves);
      
      cp.add(std::string("spiral"), std::string("ADCsPerInterleave"), 
	     (size_t)spi_parm.ADCsPerInterleave);

      cp.add(std::string("spiral"), std::string("SamplesPerADC"), 
	     (size_t)spi_parm.SamplesPerADC);

      cp.add(std::string("spiral"), std::string("SamplesToSkipStart"), 
	     (size_t)spi_parm.SamplesToSkipStart);

      cp.add(std::string("spiral"), std::string("SamplesToSkipEnd"), 
	     (size_t)spi_parm.SamplesToSkipEnd);

      cp.add(std::string("spiral"), std::string("SamplingTime_ns"), 
	     (size_t)spi_parm.SamplingTime_ns);

      cp.add(std::string("spiral"), std::string("Reordering"), 
	     (size_t)spi_parm.Reordering);
      
      cp.add(std::string("spiral"), std::string("MaxGradient_Gcm"), 
	     spi_parm.MaxGradient_Gcm);
      
      cp.add(std::string("spiral"), std::string("MaxSlewRate_Gcms"), 
	     spi_parm.MaxSlewRate_Gcms);
      
      cp.add(std::string("spiral"), std::string("krmax_cm"), 
	     spi_parm.krmax_cm);
      
      cp.add(std::string("spiral"), std::string("FOVCoeff_1"), 
	     spi_parm.FOVCoeff_1);      
    }
  }
  

  ACE_DEBUG( (LM_DEBUG, 
	      ACE_TEXT("Running with config:\n %s"), 
	      cp.exportCharArray()) );


  std::string config = cp.exportConf();

  conf.configuration_length = config.size()+1; //+1 for the null termination
  //char config_info[1024];
  


  ACE_DEBUG(( LM_INFO, ACE_TEXT("Sending configuration %s@%s \n"), 
	      conf.configurator_name, conf.configurator_lib ));
  
  ACE_INET_Addr server(port_no,hostname);
  ACE_SOCK_Connector connector;
  ACE_SOCK_Stream peer;  

  if (connector.connect(peer,server) == -1) {
    ACE_ERROR_RETURN(( LM_ERROR, ACE_TEXT("%p\n"), ACE_TEXT("connect")), 100);
  } 

  GadgetSocketReceiver mrecv(&peer);
  if (mrecv.register_reader(GADGET_MESSAGE_ACQUISITION, 
			    new GadgetAcquisitionMessageReader()) < 0) {
    ACE_ERROR_RETURN(( LM_ERROR, ACE_TEXT("%p\n"), 
		       ACE_TEXT("Unable to register acquisition reader")), -1);
  }

  if (mrecv.register_reader(GADGET_MESSAGE_IMAGE, new ImageWriter()) < 0) {
    ACE_ERROR_RETURN(( LM_ERROR, ACE_TEXT("%p\n"), 
		       ACE_TEXT("Unable to register image reader")), -1);
  }

  if (mrecv.open() < 0) {
    ACE_ERROR_RETURN(( LM_ERROR, ACE_TEXT("%p\n"), 
		       ACE_TEXT("Failed to open message receiver")), -1);
  }

  peer.send_n(&id, sizeof(GadgetMessageIdentifier));
  peer.send_n(&conf, sizeof(GadgetMessageConfigurator));
  peer.send_n(config.c_str(), conf.configuration_length);
  //peer.send_n(config_info, conf.configuration_length);

  //We need an array for collecting the data from all channels prior to transmission
  NDArray< std::complex<float> > buf;
  std::vector<int> buf_dim; 
  buf_dim.push_back(sd.GetMaxValues()->ushSamplesInScan);
  buf_dim.push_back(sd.GetMaxValues()->ushUsedChannels);

  if (!buf.create(buf_dim)) {
    ACE_DEBUG( (LM_ERROR, 
		ACE_TEXT("Unable to allocate buffer for collecting channel data\n")) );
    return -1;
  }

  //Now loop through the data and send it all to the Gadgetron
  SiemensMdhNode* next = sd.GetFirstNode();
  while (next) {
    //Now do we have any noise adjust scans
    unsigned int noise_adjust_mask = (1<<25);
    if (next->mdh.aulEvalInfoMask[0] & noise_adjust_mask) {
      std::cout << "Skipping profile" << std::endl;
      next = (SiemensMdhNode*)next->next;
      continue;
    }
 

    GadgetMessageAcquisition acq_head = acq_head_base;

    acq_head.idx.line                 = next->mdh.sLC.ushLine;
    acq_head.idx.acquisition          = next->mdh.sLC.ushAcquisition;
    acq_head.idx.slice                = next->mdh.sLC.ushSlice;
    acq_head.idx.partition            = next->mdh.sLC.ushPartition;
    acq_head.idx.echo                 = next->mdh.sLC.ushEcho;
    acq_head.idx.phase                = next->mdh.sLC.ushPhase;
    acq_head.idx.repetition           = next->mdh.sLC.ushRepetition;
    acq_head.idx.set                  = next->mdh.sLC.ushSet;
    acq_head.idx.segment              = next->mdh.sLC.ushSeg;
    acq_head.idx.channel              = next->mdh.ushChannelId;
    
    acq_head.flags |= 
      (next->mdh.aulEvalInfoMask[0] & (1 << 29)) ?
      GADGET_FLAG_LAST_ACQ_IN_SLICE : 0;
      
    acq_head.meas_uid                 = next->mdh.lMeasUID;
    acq_head.scan_counter             = next->mdh.ulScanCounter;
    acq_head.time_stamp               = next->mdh.ulTimeStamp;
    acq_head.samples                  = next->mdh.ushSamplesInScan;
    acq_head.channels                 = next->mdh.ushUsedChannels;
    acq_head.position[0]              = next->mdh.sSD.sSlicePosVec.flSag;
    acq_head.position[1]              = next->mdh.sSD.sSlicePosVec.flCor;
    acq_head.position[2]              = next->mdh.sSD.sSlicePosVec.flTra;

    memcpy(acq_head_base.quarternion, 
	   next->mdh.sSliceData.aflQuaternion,
	   sizeof(float)*4);   
    
    memcpy(buf.get_data_ptr()+next->mdh.ushChannelId*next->mdh.ushSamplesInScan,
	   next->data,
	   sizeof(float)*acq_head.samples*2);

    if (next->mdh.ushChannelId == (next->mdh.ushUsedChannels-1)) {
      //This is it, let's send it!
      id.id = GADGET_MESSAGE_ACQUISITION;

      peer.send_n(&id, sizeof(GadgetMessageIdentifier));
      peer.send_n(&acq_head, sizeof(GadgetMessageAcquisition));
      peer.send_n(buf.get_data_ptr(),
		  sizeof(float)*acq_head.samples*2*acq_head.channels );
    }

    next = (SiemensMdhNode*)next->next;
  }

  //Now we have to run the event loop and wait for images to come back
  ACE_Reactor::instance()->run_reactor_event_loop ();

  peer.close();
   
  return 0;
}
