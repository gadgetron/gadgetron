#pragma once

#include <string>
#include "core/Message.h"

namespace Gadgetron{


    class ACE_Message_Block {
      virtual ~ACE_Message_Block(){};
    };
/**
   The purpose of this case is to provide a type indepent interface to all ContainerMessages

   This interface is able to set a magic number for each type which is later on used
   instead of RTTI to "safely" cast to the right GadgetContainerMessage type

 */

template <class T> class GadgetContainerMessage : public ACE_Message_Block, private Core::TypedMessage<T>
{
  typedef ACE_Message_Block base;


public:
  /**
   *  Constructor, passing on input arguments to the contained class.
   * @param xs Variadic arguments to the contained class
   */
  template<typename... X> GadgetContainerMessage(X... xs) : Core::TypedMessage<T>(T(xs...))
   {

   }




  virtual ~GadgetContainerMessage()
  {
    //ACE_Message_Block will take care of deallocating space for the object itself;
  }

  virtual void* release()
  {    
    return nullptr; //This may be the worst line of code in the entire repository.
  }

  T* getObjectPtr() 
  {
    return this->data.get();
  }
/*
  virtual GadgetContainerMessage<T>* duplicate() 
  {
    GadgetContainerMessage<T>* nb = new GadgetContainerMessage<T>(this->data_block()->duplicate());
    nb->rd_ptr (this->rd_ptr_);
    nb->wr_ptr (this->wr_ptr_);
    if (this->cont_) {
      nb->cont_ = this->cont_->duplicate();
    }
    return nb;
  }
  */

};


template <class T> GadgetContainerMessage<T>* AsContainerMessage(ACE_Message_Block* mb)
{
  if (typeid(mb) == typeid(GadgetContainerMessage<T>*))
    return reinterpret_cast<GadgetContainerMessage<T>* >(mb);
  return nullptr;
}
}
