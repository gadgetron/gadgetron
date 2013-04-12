#ifndef GADGETCONTAINERMESSAGE_H
#define GADGETCONTAINERMESSAGE_H
#pragma once

#include <ace/Message_Block.h>
#include <string>

namespace Gadgetron{
/**
   The purpose of this case is to provide a type indepent interface to all ContainerMessages

   This interface is able to set a magic number for each type which is later on used
   instead of RTTI to "safely" cast to the right GadgetContainerMessage type

 */
class GadgetContainerMessageBase : public ACE_Message_Block
{
  typedef ACE_Message_Block base;
  
 public:

  enum { CONTAINER_MESSAGE_BLOCK = (ACE_Message_Block::USER_FLAGS << 2) };

  GadgetContainerMessageBase(size_t size) : base(size)
  {
    set_flags(CONTAINER_MESSAGE_BLOCK); //Mark this message block as a container, so that we know it is safe to type cast it.
  }

#ifdef WIN32
  std::string getTypeID() { return type_magic_id_; }
  template <class T> static std::string magic_number_for_type() { return std::string(typeid(T).name()); } 

protected:
  std::string type_magic_id_;

#else

  int getTypeID() { return type_magic_id_; }

  template <class T> static int magic_number_for_type(){
    //Will only get set once for each instanciation of this function
    static int result(next_magic_type_number()); 
    return result;
  }

 protected:
  int type_magic_id_;

  //Utility function for increting the magic number for types.
  static int next_magic_type_number()
  {
    static int magic(0);
    return magic++;
  }	 
#endif  
};

template <class T> class GadgetContainerMessage : public GadgetContainerMessageBase
{
  typedef GadgetContainerMessageBase base;

public:
  GadgetContainerMessage()
    : base(sizeof(T))
    , content_(0)
  {
    //Using placement new to put the new object at the ACE_Message_Block location
    content_ = new (this->wr_ptr()) T; 

    //Advance the write pointer appropriately.
    this->wr_ptr(sizeof(T));

    //Assign type ID that will allow us to safely cast this message.
    type_magic_id_ = magic_number_for_type<T>(); 
  }

  virtual ~GadgetContainerMessage() 
  {
    //In case the object contained in this object has allocated memory on the heap, it must be destroyed
    if (content_) content_->~T();

    //ACE_Message_Block will take care of deallocating space for the object itself;
  }

  T* getObjectPtr() 
  {
    return content_;
  }

protected:
  T* content_;
}; 

/**
   This function replaces the slower dynamic_cast which we would otherwise rely on.
   The speed of dynamic_cast varies greatly from platform to platform.

   This function is less safe since it assumes casting to ContainerMessageBase is OK
   when a certain flag is set on the ACE_Message_Block. If some user decides to use that flag
   for other purposes, it could cause major problems that are hard to debug.

   TODO: Find a more elegant solution for this.
*/
template <class T> GadgetContainerMessage<T>* AsContainerMessage(ACE_Message_Block* mb)
{
  if (!mb || !(mb->flags() & GadgetContainerMessageBase::CONTAINER_MESSAGE_BLOCK)) {
    return 0;
  }

  GadgetContainerMessageBase* mbb = reinterpret_cast<GadgetContainerMessageBase*>(mb);
  if (mbb->getTypeID() != GadgetContainerMessageBase::magic_number_for_type<T>()) {
    return 0;
  }

  return reinterpret_cast<GadgetContainerMessage<T>* >(mbb);
}
}
#endif  //GADGETCONTAINERMESSAGE_H
