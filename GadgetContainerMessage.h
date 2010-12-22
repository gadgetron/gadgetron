#ifndef GADGETCONTAINERMESSAGE_H
#define GADGETCONTAINERMESSAGE_H

#include <ace/Message_Block.h>

template <class T> class GadgetContainerMessage : public ACE_Message_Block
{
  typedef ACE_Message_Block base;

public:
  GadgetContainerMessage()
    : base(sizeof(T))
    , content_(0)
  {
    //Using placement new to put the new object at the ACE_Message_Block location
    content_ = new (this->wr_ptr()) T; 

    //Advance the write pointer appropriately.
    this->wr_ptr(sizeof(T));
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


#endif  //GADGETCONTAINERMESSAGE_H
