
#include<istream>

#include <type_traits>

namespace Gadgetron::Core::IO{


    template<class T> void read(std::istream& stream, typename std::enable_if_t<std::is_trivially_copyable_v<T>,T>& value ){
        stream.read(static_cast<char*>(&value),sizeof(value));
    }

    template<class T> void read(std::istream& stream, hoNDArray<T>& array){
        stream.read(static_cast<char*>(array.get_data_ptr()), array.get_number_of_bytes());
    }

}