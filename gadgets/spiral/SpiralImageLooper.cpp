#include "SpiralImageLooper.h"

#include "MrServers/MrVista/Ice/IceIdeaFunctors/IceSpiral/grid.h"
#include "MrServers/MrVista/Ice/IceIdeaFunctors/IceSpiral/AuxFunctions.h"
//HS: I don't have it installed
//#include "MrServers/MrImaCalcFramework/IceFramework/IceProgramSupport/FunctorUtils/IceContextExpert.h"
// for CheckAndDisplayImages
#include "MrServers/MrVista/include/Ice/IceUtils/IceUtils.h"

#include "MrServers/MrVista/Ice/IceIdeaFunctors/IceSpiral/complexFloatTools.h"

void calc_vds(double slewmax,double gradmax,double Tgsample,double Tdsample,int Ninterleaves,
			  double* fov, int numfov,double krmax,
			  int ngmax, double** xgrad,double** ygrad,int* numgrad);

void calc_traj(double* xgrad, double* ygrad, int ngrad, int Nints, double Tgsamp, double krmax,
			  double** x_trajectory, double** y_trajectory,
              double** weights);

// HS: flow
static int stflow_index=0;			  
			  




IResult SpiralImageLooper::getPhaseContrastImg(const IceAs &img0, const IceAs &img1, IceAs &imgFlow) {

	int coldim = img0.getObj()->getLen(COL);
	int lindim = img0.getObj()->getLen(LIN);
	
	int nc = img0.getObj()->getLen(CHA);
	int c,i;
	
	IceAs tmp0(img0);
	IceAs tmp1(img1);
	IceAs tmpf(imgFlow);
	
	
	for (c=0; c<nc;c++) {
		tmp0.modify(CHA, c, 1, 1);
		tmp1.modify(CHA, c, 1, 1);
		tmpf.modify(CHA, c, 1, 1);
		complexFloat* imgp0 = (complexFloat*)tmp0.calcSplObjStartAddr();
		complexFloat* imgp1 = (complexFloat*)tmp1.calcSplObjStartAddr();
		complexFloat* fp = (complexFloat*)tmpf.calcSplObjStartAddr();
		for (i=0; i<coldim*lindim;i++){			
			//fp[i] = complexFloat(angleOf(compMult(imgp0[i], conjOf(imgp1[i]))), 0.0);
			//fp[i] = compMult(imgp0[i], conjOf(imgp1[i]));
			//complexFloat phase_contrast_pixel = compMult(imgp0[i], conjOf(imgp1[i]));
			//fp[i] = complexFloat((( atan2(imagOf(phase_contrast_pixel), realOf(phase_contrast_pixel)) / M_PI + 1.0 ) / 2), 0.0) ;
			fp[i] = compMult(imgp0[i], conjOf(imgp1[i]));
		}
	}


	return I_OK;
}
	

SpiralImageLooper::SpiralImageLooper()
{
    addCB(IFirstCall);
}

SpiralImageLooper::~SpiralImageLooper()
{

}

IResult SpiralImageLooper::EndInit(IParcEnvironment* env )
{

	return I_OK;
}

IResult SpiralImageLooper::FirstCall( IceAs& srcAs, MdhProxy& aMdh, ScanControl& ctrl )
{
	//ICE_OUT("MSH: Loading data from disk");

	//Allocate some memory

	
	stflow_index=0;
	int cha_length = getNumCoils();//IceContextExpert::Instance()->getNCha();//ctrl.m_imaDimBoundaries.m_ccha;
	int slc_length = ctrl.m_imaDimBoundaries.m_cslc;
	int nsamples_per_interleave = getSamplesPerADC()*getADCsPerInterleave() 
								  - getSamplesToSkipStart() - getSamplesToSkipEnd();
	int nsamples_per_slice = nsamples_per_interleave*getNumInts();

	IceObj::Pointer co = Parc::HeapObject<IceObj>::Create();
	if (!co->create(ICE_FORMAT_FL,
					IDD, nsamples_per_slice,
					IDE, 2))
	{
		ICE_OUT("Unable to allocate trajectory in spiral image looper");
		return I_FAIL;
	}
	m_co_as = IceAs(co);
	m_co_as = (*co)();
	Ice.Preset(m_co_as,0.0);

	IceObj::Pointer we = Parc::HeapObject<IceObj>::Create();
	if (!we->create(ICE_FORMAT_FL,
					IDD, nsamples_per_slice,
					IDE, 1))
	{
		ICE_OUT("Unable to allocate trajectory in spiral image looper");
		return I_FAIL;
	}
	m_weights_as = IceAs(we);
	m_weights_as = (*we)();
	Ice.Preset(m_weights_as,0.0);

	ICE_OUT("Cha_length: ") << cha_length << std::endl;
	ICE_OUT("nsamples_per_slice: ") << nsamples_per_slice << std::endl;
	ICE_OUT("nsamples_per_interleave: ") << nsamples_per_interleave << std::endl;

	IceObj::Pointer data = Parc::HeapObject<IceObj>::Create();
	if (!data->create(ICE_FORMAT_CXFL,
					IDD, nsamples_per_slice,
					CHA, cha_length))
	{
		ICE_OUT("Unable to allocate trajectory in spiral image looper");
		return I_FAIL;
	}
	m_data_as = IceAs(data);
	m_data_as = (*data)();
	Ice.Preset(m_data_as,0.0);

	// HS: flow
	IceObj::Pointer data2 = Parc::HeapObject<IceObj>::Create();
	if (!data2->create(ICE_FORMAT_CXFL,
					IDD, nsamples_per_slice,
					CHA, cha_length))
	{
		ICE_OUT("Unable to allocate trajectory in spiral image looper");
		return I_FAIL;
	}
	m_data_as2 = IceAs(data2);
	m_data_as2 = (*data2)();
	Ice.Preset(m_data_as2,0.0);

	
	IceSDC isdc;
	isdc.calcSoda(m_co_as,srcAs);
	isdc.calcSoda(m_weights_as,srcAs);
	isdc.calcSoda(m_data_as,srcAs);
	isdc.calcSoda(m_data_as2,srcAs);
	IceSoda* sodaPtr = const_cast<IceSoda*>(m_data_as.calcCurrentSodaPtr());
	ICE_OUT("DATA SODA: ") << *sodaPtr;
#ifdef MORE_OUT
	ICE_OUT("\n\n %%%%%% Trajectory Calculation %%%%%%%%%%%%%\n\n");
#endif	
	CalculateTrajectory();
#ifdef MORE_OUT	
	ICE_OUT("\n\n %%%%%% done %%%%%%%%%%%%%\n\n");
#endif	
	//LoadDataFromDisk(srcAs);
	
	
	return I_OK;
}

IResult SpiralImageLooper::endOfJob(IResult reason)
{

	return I_OK;
}


IResult SpiralImageLooper::ComputeScan(IceAs& srcAs, MdhProxy& aMdh, ScanControl& ctrl)
{
#ifdef MORE_OUT
	ICE_OUT("\n\n <<< CheckFlow " << checkFlow(aMdh) << ">>\n\n");
#endif	
	if (getflow()) {
		ICE_OUT("\n\n<<< ComputeScanForFlowData >>>>\n\n");
		return ComputeScanForFlowData(srcAs, aMdh, ctrl);
	}
	IResult res;
    res = ExecuteCallbacks(srcAs, aMdh, ctrl);
    if(failed(res))
        ICE_RET_ERR("ExecuteCallbacks failed, aborting...", res)

	int cha_length = ctrl.m_imaDimBoundaries.m_ccha;
	int slc_length = ctrl.m_imaDimBoundaries.m_cslc;
	int nsamples_per_interleave = getSamplesPerADC()*getADCsPerInterleave() 
								  - getSamplesToSkipStart() - getSamplesToSkipEnd();
	int nsamples_per_slice = nsamples_per_interleave*getInterleaves();

	//Let's copy the incoming data to the data array
	float* src_ptr = (float*)srcAs.calcSplObjStartAddr();
	float* dst_ptr = (float*)m_data_as.calcSplObjStartAddr();
	int seg_offset;
	int samples_to_copy;
	//IceUtils::CheckAndDisplay(srcAs, IceDisplayData::DISPLAY);
	int dst_offset, src_offset;

	
	if (aMdh.isFirstScanInSlice()) {
		IceSDC isdc;
		isdc.resetSoda (m_co_as);
		isdc.resetSoda (m_weights_as);
		isdc.resetSoda (m_data_as);
		isdc.calcSoda(m_co_as,srcAs);
		isdc.calcSoda(m_weights_as,srcAs);
		isdc.calcSoda(m_data_as,srcAs);
	}

	
	for (int c = 0; c < cha_length; c++)
	{
		if (aMdh.getCseg() > 0)
		{
			samples_to_copy = getSamplesPerADC();
			seg_offset = aMdh.getCseg()*getSamplesPerADC()-getSamplesToSkipStart();
		}
		else
		{
			samples_to_copy = getSamplesPerADC()-getSamplesToSkipStart();
			seg_offset = 0;
		}

		if (aMdh.getCseg() == getADCsPerInterleave()-1)
		{
			samples_to_copy = getSamplesPerADC()-getSamplesToSkipEnd();
		}

		dst_offset = (nsamples_per_slice*c+nsamples_per_interleave*aMdh.getClin()+seg_offset)*2;
		src_offset = (c*getSamplesPerADC())*2;
		IceAs tmpAs = srcAs;
		tmpAs.modify(CHA,c,1,1);
		//ICE_OUT("srcAs") << tmpAs;
		src_ptr = (float*)tmpAs.calcSplObjStartAddr();
		/*
		ICE_OUT("nsamples_per_interleave: ") << nsamples_per_interleave;
		ICE_OUT("srcOffset: ") << (c*getSamplesPerADC())*2;
		ICE_OUT("dstOffset: ") << (nsamples_per_slice*c+nsamples_per_interleave*aMdh.getClin()+seg_offset)*2;
		ICE_OUT("Samples to copy: ") << samples_to_copy << std::endl;
		ICE_OUT("Coil: ") << c << ", " << cha_length << ", seg: " << aMdh.getCseg() << ", " << getADCsPerInterleave() << std::endl;
		*/
		for (int s = 0; s < samples_to_copy*2; s++)
		{
			dst_ptr[dst_offset+s] = src_ptr[s];
		}
	}


	if (aMdh.isLastScanInSlice() && (aMdh.getCseg() == (getADCsPerInterleave()-1)))
	{
		int dims[2];
		dims[0] = getMatrixSize();
		dims[1] = getMatrixSize();

		int fixed_dims[] = {0,0};

		//WriteData(m_co_as, "W:/temp/co_ice.real");
		//WriteData(m_weights_as, "W:/temp/weights_ice.real");
		//WriteData(m_data_as, "W:/temp/data_ice.cplx");

		IceAs gd_result;
		//LoadDataFromDisk(srcAs);
		GridData(m_co_as, m_data_as, m_weights_as, dims, fixed_dims, gd_result, true, true);

		//IceUtils::CheckAndDisplay(gd_result, IceDisplayData::DISPLAY);

		ImageControl ictrl;
		ictrl.m_imaDimBoundaries = ctrl.m_imaDimBoundaries;
		ictrl.m_imaDimBoundaries.m_ccol = getMatrixSize();
		ictrl.m_imaDimBoundaries.m_clin = getMatrixSize();
		ictrl.m_IceExtractMode = IemAmplitude;


		//Create a miniheader
		MiniHeader::Pointer imgHead = Parc::HeapObject<MiniHeader>::Create();

		IceSDC isdc;
		isdc.calcSoda(gd_result,srcAs);

		//Some scaling
		Ice.Mul(gd_result, gd_result, gd_result.getObj()->getLen(COL)*gd_result.getObj()->getLen(LIN));

		ImageReady(gd_result, imgHead, ictrl);
	}

	return I_OK;
}


IResult SpiralImageLooper::ComputeScanForFlowData(IceAs& srcAs, MdhProxy& aMdh, ScanControl& ctrl)
{

	IResult res;
    res = ExecuteCallbacks(srcAs, aMdh, ctrl);
    if(failed(res))
        ICE_RET_ERR("ExecuteCallbacks failed, aborting...", res)

	int cha_length = ctrl.m_imaDimBoundaries.m_ccha;
	int slc_length = ctrl.m_imaDimBoundaries.m_cslc;
	int nsamples_per_interleave = getSamplesPerADC()*getADCsPerInterleave() 
								  - getSamplesToSkipStart() - getSamplesToSkipEnd();
	int nsamples_per_slice = nsamples_per_interleave*getNumInts();

	//Let's copy the incoming data to the data array
	float* src_ptr = (float*)srcAs.calcSplObjStartAddr();
	float* dst_ptr = (float*)m_data_as.calcSplObjStartAddr();
	
	//HS: flow
	float* dst_ptr2 = (float*)m_data_as2.calcSplObjStartAddr();
	
	int seg_offset;
	int samples_to_copy;
	//IceUtils::CheckAndDisplay(srcAs, IceDisplayData::DISPLAY);
	int dst_offset, src_offset;

	
	if (aMdh.isFirstScanInSlice()) {
		IceSDC isdc;
		isdc.resetSoda (m_co_as);
		isdc.resetSoda (m_weights_as);
		isdc.resetSoda (m_data_as);
		isdc.resetSoda (m_data_as2);
		isdc.calcSoda(m_co_as,srcAs);
		isdc.calcSoda(m_weights_as,srcAs);
		isdc.calcSoda(m_data_as,srcAs);
		isdc.calcSoda(m_data_as2,srcAs);
	}
#ifdef MORE_OUT
	ICE_OUT("\n\n %%%%%% Copying data %%%%%%%%%%%%%\n\n");
#endif	
	for (int c = 0; c < cha_length; c++)
	{
		if (aMdh.getCseg() > 0)
		{
			samples_to_copy = getSamplesPerADC();
			seg_offset = aMdh.getCseg()*getSamplesPerADC()-getSamplesToSkipStart();
		}
		else
		{
			samples_to_copy = getSamplesPerADC()-getSamplesToSkipStart();
			seg_offset = 0;
		}

		if (aMdh.getCseg() == getADCsPerInterleave()-1)
		{
			samples_to_copy = getSamplesPerADC()-getSamplesToSkipEnd();
		}

		dst_offset = (nsamples_per_slice*c+nsamples_per_interleave*aMdh.getClin()+seg_offset)*2;
		src_offset = (c*getSamplesPerADC())*2;
		IceAs tmpAs = srcAs;
		tmpAs.modify(CHA,c,1,1);
		//ICE_OUT("srcAs") << tmpAs;
		src_ptr = (float*)tmpAs.calcSplObjStartAddr();
		/*
		ICE_OUT("nsamples_per_interleave: ") << nsamples_per_interleave;
		ICE_OUT("srcOffset: ") << (c*getSamplesPerADC())*2;
		ICE_OUT("dstOffset: ") << (nsamples_per_slice*c+nsamples_per_interleave*aMdh.getClin()+seg_offset)*2;
		ICE_OUT("Samples to copy: ") << samples_to_copy << std::endl;
		ICE_OUT("Coil: ") << c << ", " << cha_length << ", seg: " << aMdh.getCseg() << ", " << getADCsPerInterleave() << std::endl;
		*/
#if 0		
		for (int s = 0; s < samples_to_copy*2; s++)
		{
			dst_ptr[dst_offset+s] = src_ptr[s];
		}
#endif		

		// HS: flow
		int bytes_to_copy = samples_to_copy*2*sizeof(float);
		
		// find which arm first
		if (stflow_index == 0) {
#ifdef MORE_OUT		
			ICE_OUT("\n\n %%%%%% Non-Flow, interleave " << aMdh.getClin() << ", coil "<< c << "%%%%%%%%%%%%%\n\n");
#endif			
			memcpy(&dst_ptr[dst_offset], src_ptr, bytes_to_copy);
		} else { // with flow grad
#ifdef MORE_OUT		
			ICE_OUT("\n\n %%%%%% Flow, interleave " << aMdh.getClin() << ", coil "<< c << " %%%%%%%%%%%%%\n\n");
#endif			
			memcpy(&dst_ptr2[dst_offset], src_ptr, bytes_to_copy);
		}
		
			//dst_ptr[dst_offset+s] = src_ptr[s];
		
	}
#ifdef MORE_OUT	
	ICE_OUT("\n\n %%%%%% Done %%%%%%%%%%%%%\n\n");
#endif

	
	if (aMdh.isLastScanInSlice() && (aMdh.getCseg() == (getADCsPerInterleave()-1)))
	{
		
		int dims[2];
		dims[0] = getMatrixSize();
		dims[1] = getMatrixSize();

		int fixed_dims[] = {0,0};

		//WriteData(m_co_as, "W:/temp/co_ice.real");
		//WriteData(m_weights_as, "W:/temp/weights_ice.real");
		//WriteData(m_data_as, "W:/temp/data_ice.cplx");
#ifdef MORE_OUT
		ICE_OUT("\n\n %%% GRIDDING 0 %%% \n\n");
#endif		
		IceAs gd_result;
		//LoadDataFromDisk(srcAs);
		GridData(m_co_as, m_data_as, m_weights_as, dims, fixed_dims, gd_result, true, true);

#ifdef MORE_OUT		
		ICE_OUT("\n\n %%% GRIDDING 1 %%% \n\n");
#endif		
		// HS: flow
		IceAs gd_result2;
		//LoadDataFromDisk(srcAs);
		GridData(m_co_as, m_data_as2, m_weights_as, dims, fixed_dims, gd_result2, true, true);
	
		IceAs flow_result;
		int coldim = gd_result.getObj()->getLen(COL);
		int lindim = gd_result.getObj()->getLen(LIN);
		int nc = gd_result.getObj()->getLen(CHA);

#ifdef MORE_OUT		
		ICE_OUT("\n-------\n Creating object of size " << coldim << "x" << lindim<<"\n-------\n");
#endif		
		IceObj::Pointer fptr = Parc::HeapObject<IceObj>::Create();
		fptr->create(ICE_FORMAT_CXFL,
							COL, coldim,
							LIN, lindim,
							CHA, nc);
		
		flow_result = IceAs(fptr);
		flow_result = (*fptr)();	
		Ice.Preset(flow_result,0.0);
		
		//getFlowImg(gd_result, gd_result2, flow_result);
		getPhaseContrastImg(gd_result, gd_result2, flow_result);
		
		//IceUtils::CheckAndDisplay(gd_result, IceDisplayData::DISPLAY);
		

		ImageControl ictrl;
		ictrl.m_imaDimBoundaries = ctrl.m_imaDimBoundaries;
		ictrl.m_imaDimBoundaries.m_ccol = getMatrixSize();
		ictrl.m_imaDimBoundaries.m_clin = getMatrixSize();
		ictrl.m_IceExtractMode = IemAmplitude;
#ifdef MORE_OUT
		ICE_OUT("\n\n %%% Create MiniHeader %%% \n\n");
#endif		
		//Create a miniheader
		MiniHeader::Pointer imgHead = Parc::HeapObject<MiniHeader>::Create();

		//flow
		//imgHead->set("UncombinedIma", false);
#ifdef MORE_OUT		
		ICE_OUT("\n\n %%% calcSoda %%% \n\n");
#endif		
		IceSDC isdc;
		isdc.calcSoda(gd_result,srcAs);
		isdc.calcSoda(gd_result2,srcAs);
		isdc.calcSoda(flow_result,srcAs);

		ImageControl ctrl_magnitude(ctrl);
		ImageControl ctrl_neg(ctrl);
		ImageControl ctrl_pos(ctrl);
		
		//Some scaling
		//ICE_OUT("\n\n %%% Scaling %%% \n\n");
		//Ice.Mul(gd_result, gd_result, gd_result.getObj()->getLen(COL)*gd_result.getObj()->getLen(LIN));
		//Ice.Mul(flow_result, flow_result, flow_result.getObj()->getLen(COL)*flow_result.getObj()->getLen(LIN));

		//ImageReady(gd_result2, imgHead, ictrl);
		
		//IceUtils::CheckAndDisplay(flow_result, IceDisplayData::DISPLAY);
		ImageReady(flow_result, imgHead, ictrl);
		
	} 
	
	stflow_index = !stflow_index;
	return I_OK;
}

IResult SpiralImageLooper::CalculateTrajectory()
{
    double  smax = getMaxSlewRate_Gcms();   	/*	Maximum slew rate, G/cm/s		*/
    double  gmax = getMaxGradient_Gcm();       /* 	maximum gradient amplitude, G/cm	*/
    double  Tdsamp = getSamplingTime_ns()/(1.0e9);     /*	Data Sample period (s)			*/
    double  Tgsamp = 1e-5;     /*	Data Sample period (s)			*/
    int     Nints = getNumInts();  	/*	Number of interleaves			*/
    double  fov = getFOVCoeff_1();	/*	FOV coefficients		*/
    int     Nfov = 1;       /*	Number of FOV coefficients		*/
    double  krmax = getkrmax_cm();  	/*	Maximum k-space extent (/cm)		*/
    int     ngmax = 1e5;      /*	Maximum number of gradient samples	*/
    double  *xgrad; 	/* 	X-component of gradient.	*/
    double  *ygrad;     /*	Y-component of gradient.	*/
    double  *x_trajectory;
    double  *y_trajectory;
    double  *weighting;
    int     ngrad;	
    int     count;

#ifdef MORE_OUT	
	ICE_OUT("Calculating Trajectory: ") << std::endl 
		<< "S_max: " << smax << std::endl
		<< "G_max: " << gmax << std::endl
		<< "Tdsamp: " << Tdsamp << std::endl
		<< "Nints:  " << Nints << std::endl
		<< "fov: " << fov << std::endl
		<< "krmax: " << krmax << std::endl;
#endif
    /*	Call C-function Here to calculate gradients */
    calc_vds(smax,gmax,Tdsamp,Tdsamp,Nints,&fov,Nfov,krmax,ngmax,
            &xgrad,&ygrad,&ngrad);
 
	ICE_OUT("MSH: ngrad: ") << ngrad;
	if ((getADCsPerInterleave()*getSamplesPerADC()-ngrad) != getSamplesToSkipEnd())
	{
		ngrad = (getADCsPerInterleave()*getSamplesPerADC()-getSamplesToSkipEnd()) < ngrad ? getADCsPerInterleave()*getSamplesPerADC()-getSamplesToSkipEnd() : ngrad;
		setSamplesToSkipEnd(getADCsPerInterleave()*getSamplesPerADC()-ngrad);
	}
	ICE_OUT("MSH: ngrad: ") << ngrad;

    /* Calcualte the trajectory and weights*/
    calc_traj(xgrad, ygrad, ngrad, Nints, Tdsamp, krmax, &x_trajectory, &y_trajectory, &weighting);

	float* co_ptr = (float*)m_co_as.calcSplObjStartAddr();
	float* we_ptr = (float*)m_weights_as.calcSplObjStartAddr();

	float xmax = 0;
    for (int i = 0; i < (ngrad*Nints); i++)
    {
		if (x_trajectory[i] > xmax) xmax = x_trajectory[i];
        co_ptr[i] = -x_trajectory[i]*(getMatrixSize()>>1);
        co_ptr[i+(ngrad*Nints)] = -y_trajectory[i]*(getMatrixSize()>>1);
        we_ptr[i] = weighting[i];
    }
	ICE_OUT("xmax: ") << xmax << std::endl;

    delete [] xgrad;
    delete [] ygrad;
    delete [] x_trajectory;
    delete [] y_trajectory;
    delete [] weighting;

	return I_OK;
}

IResult SpiralImageLooper::LoadDataFromDisk(IceAs& srcAs)
{
	int ndim;
	int dims[10];
	FILE* f;

	//Load trajectory
	/*
	f = fopen("c:/MSH/spiral/co.real","rb");

	fread(&ndim,sizeof(int),1,f);
	fread(&dims,sizeof(int),ndim,f);
	fread((float*)m_co_as.calcSplObjStartAddr(),sizeof(float),dims[0]*dims[1],f);
	fclose(f);
	*/
	//Load weights
	f = fopen("c:/MSH/spiral/weights.real","rb");

	fread(&ndim,sizeof(int),1,f);
	fread(&dims,sizeof(int),ndim,f);
	fread((float*)m_weights_as.calcSplObjStartAddr(),sizeof(float),dims[0],f);
	fclose(f);

	//Load data
	/*
	f = fopen("c:/MSH/spiral/data.cplx","rb");

	fread(&ndim,sizeof(int),1,f);
	fread(&dims,sizeof(int),ndim,f);
	fread((float*)m_data_as.calcSplObjStartAddr(),sizeof(float),2*dims[0],f);
	fclose(f);
	*/

	ICE_OUT("Data loaded from disk:");
	ICE_OUT("m_co_as: ") << m_co_as;
	ICE_OUT("m_weights_as: ") << m_weights_as;
	ICE_OUT("m_data_as: ") << m_data_as;


	return I_OK;
}




IResult SpiralImageLooper::GridData(IceAs& trajectory, IceAs& data, IceAs& weights, int* dims, int* fixed_dims, IceAs& result, bool result_in_imagespace, bool remove_padding)
{
	long ndim = trajectory.getObj()->getLen(IDE);
	long cha_length = data.getObj()->getLen(CHA);

	gridding_parameters parm;
	point_vectors pv;
	
	InitGridding(trajectory, dims, fixed_dims, parm, pv);

	ScaleGriddingWeights(weights, parm, ndim, dims, fixed_dims);

	int oversampled_dimensions[3];
	long points = 1;
	for (int i = 0; i < ndim; i++)
	{
		oversampled_dimensions[i] = (int)(dims[i]+(parm.over_sampling-1.0)*dims[i]*(!fixed_dims[i]));
		points *= oversampled_dimensions[i];
	}


	IceAs deap;
	CalculateDeapodization(parm,pv,dims,fixed_dims,deap);
	Ice.Mul(deap,deap,points); 

	IceAs gd_r;
	GridConvolution(parm,pv,data,weights,dims,fixed_dims,1,gd_r);

	K2I(gd_r,COL);
	K2I(gd_r,LIN);
	if (ndim > 2) K2I(gd_r,PHS); //Should be more general in the fufure

	ApplyDeapodization(gd_r, deap);

 	IceObj::Pointer res = Parc::HeapObject<IceObj>::Create();

	if (remove_padding)
	{
		IceAs copy_as(gd_r);
		if (ndim == 2)
		{
			res->create(ICE_FORMAT_CXFL,
						COL, dims[0],
						LIN, dims[1],
						CHA, cha_length);

			copy_as.init(COL,(oversampled_dimensions[0]-dims[0])>>1, dims[0],1,
						 LIN,(oversampled_dimensions[1]-dims[1])>>1, dims[1],1,
						 CHA,0,cha_length,1);	
		}
		else if (ndim == 3)
		{
			res->create(ICE_FORMAT_CXFL,
						COL, dims[0],
						LIN, dims[1],
						PHS, dims[2],
						CHA, cha_length); //PHS, This is hard coded for now, should be flexible....TODO

			copy_as.init(COL,(oversampled_dimensions[0]-dims[0])>>1, dims[0],1,
						 LIN,(oversampled_dimensions[1]-dims[1])>>1, dims[1],1,
						 PHS,(oversampled_dimensions[2]-dims[2])>>1, dims[2],1,
						 CHA,0,cha_length,1);
		}

		result = IceAs(res);
		result = (*res)();

		Ice.Copy(result,copy_as);
	}
	else
	{
		result = gd_r;
	}

	if (!result_in_imagespace)
	{
		if (ndim > 2) I2K(result,PHS); //Should be more general in the fufure

		I2K(result,COL);
		I2K(result,LIN);
	}

	return I_OK;
}

IResult SpiralImageLooper::GridConvolution(gridding_parameters& parm, point_vectors& pv, IceAs& data, IceAs& weights, int* dims, int* fixed_dims, int direction, IceAs& result)
{	
	int npoints    = pv.number_of_points;
	int dimensions = pv.number_of_dimensions;
	int channels   = data.getObj()->getLen(CHA);

	IceObj::Pointer res = Parc::HeapObject<IceObj>::Create();

	if (direction == 1)
	{
		if (dimensions == 2)
		{
			res->create(ICE_FORMAT_CXFL,
						COL, (int)(dims[0]+(parm.over_sampling-1.0)*dims[0]*(!fixed_dims[0])),
						LIN, (int)(dims[1]+(parm.over_sampling-1.0)*dims[1]*(!fixed_dims[1])),
						CHA, channels);
		}
		else if (dimensions == 3)
		{
			res->create(ICE_FORMAT_CXFL,
						COL, (int)(dims[0]+(parm.over_sampling-1.0)*dims[0]*(!fixed_dims[0])),
						LIN, (int)(dims[1]+(parm.over_sampling-1.0)*dims[1]*(!fixed_dims[1])),
						PHS, (int)(dims[2]+(parm.over_sampling-1.0)*dims[2]*(!fixed_dims[2])), //This is hard coded for now, should be flexible....TODO
						CHA, channels); 
		}
		else
		{
			ICE_OUT("Gridding only work in max 3 dimensions for now");
			return I_FAIL;
		}
	} 
	else
	{
		res->create(ICE_FORMAT_CXFL,
					IDD, npoints,
					CHA, channels);
	}

	result = IceAs(res);
	result = (*res)();
	
	Ice.Preset(result,0.0);

	IceSDC isdc;
	isdc.calcSoda(result,data);

	float* weights_adr	= NULL;
	if (direction == 1) weights_adr = (float*)weights.calcSplObjStartAddr();

	int cha;
	IceAs tmp1(data);
	IceAs tmp2(result);
	for (cha = 0; cha < channels; cha++)
	{
		tmp1 = (*data.getObj())(CHA,cha);
		tmp2 = (*result.getObj())(CHA,cha);

		float* data_in_adr	= (float*)tmp1.calcSplObjStartAddr();
		float* data_out_adr = (float*)tmp2.calcSplObjStartAddr();
		if (dimensions == 2)
		{
			genkt_grid_2d(&parm,pv.kernel_positions,pv.grid_positions,npoints,data_in_adr,weights_adr,data_out_adr,dims,direction,fixed_dims);
		}
		else if (dimensions == 3)
		{
			genkt_grid_3d(&parm,pv.kernel_positions,pv.grid_positions,npoints,data_in_adr,weights_adr,data_out_adr,dims,direction,fixed_dims);
		}
		else
		{
			ICE_OUT("Gridding only works in max 3 dimensions for now");
			return I_FAIL;
		}
	}

	return I_OK;
}

IResult SpiralImageLooper::InitGridding(IceAs& trajectory, int* dims, int* fixed_dims, gridding_parameters& parm, point_vectors& pv)
{

	int npoints    = trajectory.getObj()->getLen(IDD);
	int dimensions = trajectory.getObj()->getLen(IDE);

	init_gridding_parameters(&parm);

	calculate_kernel_tables(&parm);
	calculate_point_vectors(trajectory, dims, &parm, &pv, fixed_dims);

	return I_OK;
}


IResult SpiralImageLooper::CalculateDeapodization(gridding_parameters& parm, point_vectors& pv, int* dims, int* fixed_dims, IceAs& result)
{
	
	int ndim = pv.number_of_dimensions;

	int oversampled_dimensions[3];
	for (int i = 0; i < ndim; i++)
	{
		oversampled_dimensions[i] = (int)(dims[i]+(parm.over_sampling-1.0)*dims[i]*(!fixed_dims[i]));
	}

  	IceObj::Pointer res = Parc::HeapObject<IceObj>::Create();
	if (ndim == 2)
	{
		res->create(ICE_FORMAT_CXFL,
					COL, oversampled_dimensions[0],
					LIN, oversampled_dimensions[1]);
	}
	else if (ndim == 3)
	{
		res->create(ICE_FORMAT_CXFL,
					COL, oversampled_dimensions[0],
					LIN, oversampled_dimensions[1],
					PHS, oversampled_dimensions[2]); //PHS, This is hard coded for now, should be flexible....TODO
	}
	else
	{
		ICE_OUT("Deapodization: Gridding only work in max 3 dimensions for now");
		return I_FAIL;
	}
	result = IceAs(res);
	result = (*res)();
	
	Ice.Preset(result,0.0);

	IceSDC isdc;
	isdc.calcSoda(result,pv.kernel_positions_as);

	int kernel_table_center = ((static_cast<int>(parm.kernel_table_steps*parm.kernel_width))>>1);
  
	int zmin, zmax, ymin, ymax, xmin, xmax, zoffset, yoffset, xoffset;

	if (ndim == 3 && fixed_dims[2] == 0)
	{
		zmin = static_cast<int>(-floor(parm.kernel_width/2)*parm.over_sampling);
		zmax = static_cast<int>(floor(parm.kernel_width/2)*parm.over_sampling);
	}
	else
	{
		zmin = 0; zmax = 0;
	}
	if (ndim == 3)
	{
		zoffset = oversampled_dimensions[2]>>1;
	}
	else
	{
		zoffset = 0;
	}

	if (fixed_dims[1] == 0)
	{
		ymin = static_cast<int>(-floor(parm.kernel_width/2)*parm.over_sampling);
		ymax = static_cast<int>(floor(parm.kernel_width/2)*parm.over_sampling);
	}
	else
	{
		ymin = 0; ymax = 0;
	}
	yoffset = oversampled_dimensions[1]>>1;

	if (fixed_dims[0] == 0)
	{
		xmin = static_cast<int>(-floor(parm.kernel_width/2)*parm.over_sampling);
		xmax = static_cast<int>(floor(parm.kernel_width/2)*parm.over_sampling);
	}
	else
	{
		xmin = 0; xmax = 0;
	}
	xoffset = oversampled_dimensions[0]>>1;


	double kz, ky, kx;
	int points = 0;
	int ix,iy,iz;
	float* a = (float*)result.calcSplObjStartAddr();

	for (int z = zmin; z <= zmax; z++)
	{
		iz = static_cast<int>(kernel_table_center+(z*(parm.kernel_table_steps/parm.over_sampling)));
		if (ndim == 3 && !fixed_dims[2])
		{
		  kz = parm.kernel_table[iz];
		}
		else
		{
		  kz = 1;
		}
		for (int y = ymin; y <= ymax; y++)
		{
			iy = static_cast<int>(kernel_table_center+(y*(parm.kernel_table_steps/parm.over_sampling)));
			ky = parm.kernel_table[iy];
			for (int x = xmin; x <= xmax; x++)
			{
				ix = static_cast<int>(kernel_table_center+(x*(parm.kernel_table_steps/parm.over_sampling)));
				kx = parm.kernel_table[ix];
				a[((z+zoffset)*oversampled_dimensions[1]*oversampled_dimensions[0]+(y+yoffset)*oversampled_dimensions[0]+(x+xoffset))*2] = (float)(kx*ky*kz);
				points++;
			}
		}
	}

	isdc.deactivateNextSDC();
	K2I(result,COL);
	isdc.deactivateNextSDC();
	K2I(result,LIN);
	if (ndim == 3)
	{
		isdc.deactivateNextSDC();
		K2I(result,PHS);
	}

	return I_OK;
}

IResult SpiralImageLooper::ScaleGriddingWeights(IceAs& weights, gridding_parameters& parm, int ndim, int* dims, int* fixed_dims)
{
	long npoints = weights.getObj()->getLen(IDD);

	float* w_adr = (float*)weights.calcSplObjStartAddr();

	double scale = 0;
	long i;
	for (i = 0; i < npoints; i++)
	{
		scale += w_adr[i];
	}

	int oversampled_dimensions[3];
	for (i = 0; i < ndim; i++)
	{
		oversampled_dimensions[i] = (int)(dims[i]+(parm.over_sampling-1.0)*dims[i]*(!fixed_dims[i]));
	}

	double area = oversampled_dimensions[0]*oversampled_dimensions[1]*(M_PI/4);

	if (ndim == 3)
	{
		area *= oversampled_dimensions[2];
	}

	scale /= area;
	for (i = 0; i < npoints; i++)
	{
		w_adr[i] /= (float)scale;
	}

	return I_OK;
}


IResult SpiralImageLooper::ApplyDeapodization(IceAs& data, IceAs& deap, bool inverse)
{

	int cha_length = data.getObj()->getLen(CHA) > 1 ? data.getObj()->getLen(CHA) : 1;
	IceAs iterAs(data);

	int cha;
	for (cha = 0; cha < cha_length; cha++)
	{
		iterAs = (*data.getObj())(CHA,cha);
//		if (inverse)
//		{
//			Ice.Mul(iterAs,iterAs,deap);
//		}
//		else
//		{
			Ice.Div(iterAs,iterAs,deap);
//		}
	}

	return I_OK;
}
