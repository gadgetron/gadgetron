
namespace Gadgetron
{

template <typename T> 
hoNDArrayMemoryManaged<T>::hoNDArrayMemoryManaged()
{
    BaseClass();
}

template <typename T> 
hoNDArrayMemoryManaged<T>::hoNDArrayMemoryManaged(MemManagerType& mem_manager)
: mem_manager_(mem_manager)
{
    BaseClass();
}

template <typename T> 
hoNDArrayMemoryManaged<T>::hoNDArrayMemoryManaged(std::vector<unsigned long long> *dimensions, MemManagerType& mem_manager)
: mem_manager_(mem_manager)
{
    this->create(dimensions);
}

template <typename T> 
hoNDArrayMemoryManaged<T>::hoNDArrayMemoryManaged(unsigned long long len, MemManagerType& mem_manager)
: mem_manager_(mem_manager)
{
    std::vector<unsigned long long> dim(1);
    dim[0] = len;
    this->create(dim);
}

template <typename T> 
hoNDArrayMemoryManaged<T>::hoNDArrayMemoryManaged(unsigned long long sx, unsigned long long sy, MemManagerType& mem_manager)
: mem_manager_(mem_manager)
{
    std::vector<unsigned long long> dim(2);
    dim[0] = sx;
    dim[1] = sy;
    this->create(dim);
}

template <typename T> 
hoNDArrayMemoryManaged<T>::hoNDArrayMemoryManaged(unsigned long long sx, unsigned long long sy, unsigned long long sz, MemManagerType& mem_manager)
: mem_manager_(mem_manager)
{
    std::vector<unsigned long long> dim(3);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    this->create(dim);
}

template <typename T> 
hoNDArrayMemoryManaged<T>::hoNDArrayMemoryManaged(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, MemManagerType& mem_manager)
: mem_manager_(mem_manager)
{
    std::vector<unsigned long long> dim(4);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    this->create(dim);
}

template <typename T> 
hoNDArrayMemoryManaged<T>::hoNDArrayMemoryManaged(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, MemManagerType& mem_manager)
: mem_manager_(mem_manager)
{
    std::vector<unsigned long long> dim(5);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    dim[4] = sp;
    this->create(dim);
}

template <typename T> 
hoNDArrayMemoryManaged<T>::hoNDArrayMemoryManaged(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, MemManagerType& mem_manager)
: mem_manager_(mem_manager)
{
    std::vector<unsigned long long> dim(6);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    dim[4] = sp;
    dim[5] = sq;
    this->create(dim);
}

template <typename T> 
hoNDArrayMemoryManaged<T>::hoNDArrayMemoryManaged(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, MemManagerType& mem_manager)
: mem_manager_(mem_manager)
{
    std::vector<unsigned long long> dim(7);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    dim[4] = sp;
    dim[5] = sq;
    dim[6] = sr;
    this->create(dim);
}

template <typename T> 
hoNDArrayMemoryManaged<T>::hoNDArrayMemoryManaged(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, unsigned long long ss, MemManagerType& mem_manager)
: mem_manager_(mem_manager)
{
    std::vector<unsigned long long> dim(8);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    dim[4] = sp;
    dim[5] = sq;
    dim[6] = sr;
    dim[7] = ss;
    this->create(dim);
}

template <typename T> 
hoNDArrayMemoryManaged<T>::hoNDArrayMemoryManaged(boost::shared_ptr< std::vector<unsigned long long> > dimensions, MemManagerType& mem_manager)
: mem_manager_(mem_manager)
{
    this->create(dimensions.get());
}

template <typename T> 
hoNDArrayMemoryManaged<T>::~hoNDArrayMemoryManaged()
{
    // GADGET_MSG("In ~hoNDArrayMemoryManaged()");
    if (this->delete_data_on_destruct_)
    {
        //if ( mem_manager_ )
        //{
            this->deallocate_memory();
            this->data_ = NULL;
        //}
    }
}

template <typename T> 
hoNDArrayMemoryManaged<T>::hoNDArrayMemoryManaged(const hoNDArrayMemoryManaged<T>& a)
{
    this->mem_manager_ = a.mem_manager_;
    this->data_ = 0;
    this->dimensions_ = a.dimensions_;
    this->offsetFactors_ = a.offsetFactors_;
    this->allocate_memory();
    memcpy( this->data_, a.data_, this->elements_*sizeof(T) );
}

template <typename T> 
hoNDArrayMemoryManaged<T>::hoNDArrayMemoryManaged(const hoNDArray<T>& a, MemManagerType& mem_manager)
: mem_manager_(mem_manager)
{
    this->data_ = 0;
    this->dimensions_ = a.get_dimensions();
    this->offsetFactors_ = a.get_offset_factor();
    this->allocate_memory();
    memcpy( this->data_, const_cast<T*>(a.begin()), this->elements_*sizeof(T) );
}

template <typename T> 
hoNDArrayMemoryManaged<T>& hoNDArrayMemoryManaged<T>::operator=(const hoNDArrayMemoryManaged<T>& rhs)
{
    if ( &rhs == this ) return *this;
    BaseClass::operator=( dynamic_cast< const hoNDArray<T>& >(rhs) );
    return *this;
}

template <typename T> 
hoNDArrayMemoryManaged<T>& hoNDArrayMemoryManaged<T>::operator=(const hoNDArray<T>& rhs)
{
    if ( &rhs == this ) return *this;
    BaseClass::operator=( rhs );
    return *this;
}

template <typename T> 
void hoNDArrayMemoryManaged<T>::setMemoryManager(MemManagerType& mem_manager)
{
    mem_manager_ = mem_manager;
}

template <typename T> 
void hoNDArrayMemoryManaged<T>::print(std::ostream& os) const
{
    using namespace std;

    os.unsetf(std::ios::scientific);
    os.setf(ios::fixed);

    unsigned long long i;

    os << "-------------- Gagdgetron ND Array controlled by GtPlus memory manager -------------" << endl;
    os << "Array dimension is : " << dimensions_->size() << endl;

    os << "Array size is : ";
    for (i=0; i<dimensions_->size(); i++ ) 
        os << (*dimensions_)[i] << " "; 
    os << endl;

    int elemTypeSize = sizeof(T);
    std::string elemTypeName = std::string(typeid(T).name());

    os << "Array data type is : " << elemTypeName << std::endl;
    os << "Byte number for each element is : " << elemTypeSize << std::endl;
    os << "Number of array size in bytes is : ";
    os << elements_*elemTypeSize << std::endl;

    os << std::endl;
}

template <typename T> 
inline void hoNDArrayMemoryManaged<T>::_allocate_memory( unsigned long long size, float** data )
{
    if ( mem_manager_ )
    {
        *data = (float*) mem_manager_->allocate(size*sizeof(float));
    }
    else
    {
        BaseClass::_allocate_memory(size, data);
    }
}

template <typename T> 
inline void hoNDArrayMemoryManaged<T>::_deallocate_memory( float* data )
{
    if ( mem_manager_ )
    {
        mem_manager_->free( (void*)data );
    }
    else
    {
        BaseClass::_deallocate_memory(data);
    }
}

template <typename T> 
inline void hoNDArrayMemoryManaged<T>::_allocate_memory( unsigned long long size, double** data )
{
    if ( mem_manager_ )
    {
        *data = (double*) mem_manager_->allocate(size*sizeof(double));
    }
    else
    {
        BaseClass::_allocate_memory(size, data);
    }
}

template <typename T> 
inline void hoNDArrayMemoryManaged<T>::_deallocate_memory( double* data )
{
    if ( mem_manager_ )
    {
        mem_manager_->free( (void*)data );
    }
    else
    {
        BaseClass::_deallocate_memory(data);
    }
}

template <typename T> 
inline void hoNDArrayMemoryManaged<T>::_allocate_memory( unsigned long long size, std::complex<float>** data )
{
    //GADGET_MSG("-----> In hoNDArrayMemoryManaged::_allocate_memory(size) : " << size/1024/1024.0 << " MegaBytes ");
    if ( mem_manager_ )
    {
        *data = (std::complex<float>*) mem_manager_->allocate(size*sizeof(std::complex<float>));
    }
    else
    {
        BaseClass::_allocate_memory(size, data);
    }
}

template <typename T> 
inline void hoNDArrayMemoryManaged<T>::_deallocate_memory( std::complex<float>* data )
{
    if ( mem_manager_ )
    {
        mem_manager_->free( (void*)data );
    }
    else
    {
        BaseClass::_deallocate_memory(data);
    }
}

template <typename T> 
inline void hoNDArrayMemoryManaged<T>::_allocate_memory( unsigned long long size, std::complex<double>** data )
{
    if ( mem_manager_ )
    {
        *data = (std::complex<double>*) mem_manager_->allocate(size*sizeof(std::complex<double>));
    }
    else
    {
        BaseClass::_allocate_memory(size, data);
    }
}

template <typename T> 
inline void hoNDArrayMemoryManaged<T>::_deallocate_memory( std::complex<double>* data )
{
    if ( mem_manager_ )
    {
        mem_manager_->free( (void*)data );
    }
    else
    {
        BaseClass::_deallocate_memory(data);
    }
}

template <typename T> 
inline void hoNDArrayMemoryManaged<T>::_allocate_memory( unsigned long long size, float_complext** data )
{
    if ( mem_manager_ )
    {
        *data = (float_complext*) mem_manager_->allocate(size*sizeof(float_complext));
    }
    else
    {
        BaseClass::_allocate_memory(size, data);
    }
}

template <typename T> 
inline void hoNDArrayMemoryManaged<T>::_deallocate_memory( float_complext* data )
{
    if ( mem_manager_ )
    {
        mem_manager_->free( (void*)data );
    }
    else
    {
        BaseClass::_deallocate_memory(data);
    }
}

template <typename T> 
inline void hoNDArrayMemoryManaged<T>::_allocate_memory( unsigned long long size, double_complext** data )
{
    if ( mem_manager_ )
    {
        *data = (double_complext*) mem_manager_->allocate(size*sizeof(double_complext));
    }
    else
    {
        BaseClass::_allocate_memory(size, data);
    }
}

template <typename T> 
inline void hoNDArrayMemoryManaged<T>::_deallocate_memory( double_complext* data )
{
    if ( mem_manager_ )
    {
        mem_manager_->free( (void*)data );
    }
    else
    {
        BaseClass::_deallocate_memory(data);
    }
}

}

