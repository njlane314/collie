// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME iodIsrcdIdIcollieio_dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "io/include//CollieDistribution.hh"
#include "io/include//CollieEventList.hh"
#include "io/include//CollieChannel.hh"
#include "io/include//CollieIterator.hh"
#include "io/include//CollieHistogram.hh"
#include "io/include//CollieQuickKEYS.hh"
#include "io/include//CollieIOFile.hh"
#include "io/include//CollieXsec.hh"
#include "io/include//CollieMasspoint.hh"
#include "io/include//CDF_DZero_Distribution.hh"
#include "io/include//CDF_DZero_IOpoint.hh"
#include "io/include//CDF_DZero_IOfile.hh"
#include "io/include//InfoForEfficiencyCalculation.hh"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *CollieHistogram_Dictionary();
   static void CollieHistogram_TClassManip(TClass*);
   static void *new_CollieHistogram(void *p = 0);
   static void *newArray_CollieHistogram(Long_t size, void *p);
   static void delete_CollieHistogram(void *p);
   static void deleteArray_CollieHistogram(void *p);
   static void destruct_CollieHistogram(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CollieHistogram*)
   {
      ::CollieHistogram *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::CollieHistogram));
      static ::ROOT::TGenericClassInfo 
         instance("CollieHistogram", "CollieHistogram.hh", 71,
                  typeid(::CollieHistogram), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &CollieHistogram_Dictionary, isa_proxy, 0,
                  sizeof(::CollieHistogram) );
      instance.SetNew(&new_CollieHistogram);
      instance.SetNewArray(&newArray_CollieHistogram);
      instance.SetDelete(&delete_CollieHistogram);
      instance.SetDeleteArray(&deleteArray_CollieHistogram);
      instance.SetDestructor(&destruct_CollieHistogram);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CollieHistogram*)
   {
      return GenerateInitInstanceLocal((::CollieHistogram*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::CollieHistogram*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *CollieHistogram_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::CollieHistogram*)0x0)->GetClass();
      CollieHistogram_TClassManip(theClass);
   return theClass;
   }

   static void CollieHistogram_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_InfoForEfficiencyCalculation(void *p = 0);
   static void *newArray_InfoForEfficiencyCalculation(Long_t size, void *p);
   static void delete_InfoForEfficiencyCalculation(void *p);
   static void deleteArray_InfoForEfficiencyCalculation(void *p);
   static void destruct_InfoForEfficiencyCalculation(void *p);
   static void streamer_InfoForEfficiencyCalculation(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::InfoForEfficiencyCalculation*)
   {
      ::InfoForEfficiencyCalculation *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::InfoForEfficiencyCalculation >(0);
      static ::ROOT::TGenericClassInfo 
         instance("InfoForEfficiencyCalculation", ::InfoForEfficiencyCalculation::Class_Version(), "InfoForEfficiencyCalculation.hh", 28,
                  typeid(::InfoForEfficiencyCalculation), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::InfoForEfficiencyCalculation::Dictionary, isa_proxy, 16,
                  sizeof(::InfoForEfficiencyCalculation) );
      instance.SetNew(&new_InfoForEfficiencyCalculation);
      instance.SetNewArray(&newArray_InfoForEfficiencyCalculation);
      instance.SetDelete(&delete_InfoForEfficiencyCalculation);
      instance.SetDeleteArray(&deleteArray_InfoForEfficiencyCalculation);
      instance.SetDestructor(&destruct_InfoForEfficiencyCalculation);
      instance.SetStreamerFunc(&streamer_InfoForEfficiencyCalculation);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::InfoForEfficiencyCalculation*)
   {
      return GenerateInitInstanceLocal((::InfoForEfficiencyCalculation*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::InfoForEfficiencyCalculation*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_CollieDistribution(void *p = 0);
   static void *newArray_CollieDistribution(Long_t size, void *p);
   static void delete_CollieDistribution(void *p);
   static void deleteArray_CollieDistribution(void *p);
   static void destruct_CollieDistribution(void *p);
   static void streamer_CollieDistribution(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CollieDistribution*)
   {
      ::CollieDistribution *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::CollieDistribution >(0);
      static ::ROOT::TGenericClassInfo 
         instance("CollieDistribution", ::CollieDistribution::Class_Version(), "CollieDistribution.hh", 36,
                  typeid(::CollieDistribution), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::CollieDistribution::Dictionary, isa_proxy, 16,
                  sizeof(::CollieDistribution) );
      instance.SetNew(&new_CollieDistribution);
      instance.SetNewArray(&newArray_CollieDistribution);
      instance.SetDelete(&delete_CollieDistribution);
      instance.SetDeleteArray(&deleteArray_CollieDistribution);
      instance.SetDestructor(&destruct_CollieDistribution);
      instance.SetStreamerFunc(&streamer_CollieDistribution);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CollieDistribution*)
   {
      return GenerateInitInstanceLocal((::CollieDistribution*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::CollieDistribution*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_CollieEventList(void *p = 0);
   static void *newArray_CollieEventList(Long_t size, void *p);
   static void delete_CollieEventList(void *p);
   static void deleteArray_CollieEventList(void *p);
   static void destruct_CollieEventList(void *p);
   static void streamer_CollieEventList(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CollieEventList*)
   {
      ::CollieEventList *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::CollieEventList >(0);
      static ::ROOT::TGenericClassInfo 
         instance("CollieEventList", ::CollieEventList::Class_Version(), "CollieEventList.hh", 25,
                  typeid(::CollieEventList), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::CollieEventList::Dictionary, isa_proxy, 16,
                  sizeof(::CollieEventList) );
      instance.SetNew(&new_CollieEventList);
      instance.SetNewArray(&newArray_CollieEventList);
      instance.SetDelete(&delete_CollieEventList);
      instance.SetDeleteArray(&deleteArray_CollieEventList);
      instance.SetDestructor(&destruct_CollieEventList);
      instance.SetStreamerFunc(&streamer_CollieEventList);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CollieEventList*)
   {
      return GenerateInitInstanceLocal((::CollieEventList*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::CollieEventList*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_CollieMasspoint(void *p = 0);
   static void *newArray_CollieMasspoint(Long_t size, void *p);
   static void delete_CollieMasspoint(void *p);
   static void deleteArray_CollieMasspoint(void *p);
   static void destruct_CollieMasspoint(void *p);
   static void streamer_CollieMasspoint(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CollieMasspoint*)
   {
      ::CollieMasspoint *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::CollieMasspoint >(0);
      static ::ROOT::TGenericClassInfo 
         instance("CollieMasspoint", ::CollieMasspoint::Class_Version(), "CollieMasspoint.hh", 22,
                  typeid(::CollieMasspoint), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::CollieMasspoint::Dictionary, isa_proxy, 16,
                  sizeof(::CollieMasspoint) );
      instance.SetNew(&new_CollieMasspoint);
      instance.SetNewArray(&newArray_CollieMasspoint);
      instance.SetDelete(&delete_CollieMasspoint);
      instance.SetDeleteArray(&deleteArray_CollieMasspoint);
      instance.SetDestructor(&destruct_CollieMasspoint);
      instance.SetStreamerFunc(&streamer_CollieMasspoint);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CollieMasspoint*)
   {
      return GenerateInitInstanceLocal((::CollieMasspoint*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::CollieMasspoint*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_CollieChannel(void *p = 0);
   static void *newArray_CollieChannel(Long_t size, void *p);
   static void delete_CollieChannel(void *p);
   static void deleteArray_CollieChannel(void *p);
   static void destruct_CollieChannel(void *p);
   static void streamer_CollieChannel(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CollieChannel*)
   {
      ::CollieChannel *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::CollieChannel >(0);
      static ::ROOT::TGenericClassInfo 
         instance("CollieChannel", ::CollieChannel::Class_Version(), "CollieChannel.hh", 25,
                  typeid(::CollieChannel), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::CollieChannel::Dictionary, isa_proxy, 16,
                  sizeof(::CollieChannel) );
      instance.SetNew(&new_CollieChannel);
      instance.SetNewArray(&newArray_CollieChannel);
      instance.SetDelete(&delete_CollieChannel);
      instance.SetDeleteArray(&deleteArray_CollieChannel);
      instance.SetDestructor(&destruct_CollieChannel);
      instance.SetStreamerFunc(&streamer_CollieChannel);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CollieChannel*)
   {
      return GenerateInitInstanceLocal((::CollieChannel*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::CollieChannel*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *CollieIterator_Dictionary();
   static void CollieIterator_TClassManip(TClass*);
   static void delete_CollieIterator(void *p);
   static void deleteArray_CollieIterator(void *p);
   static void destruct_CollieIterator(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CollieIterator*)
   {
      ::CollieIterator *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::CollieIterator));
      static ::ROOT::TGenericClassInfo 
         instance("CollieIterator", "CollieIterator.hh", 10,
                  typeid(::CollieIterator), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &CollieIterator_Dictionary, isa_proxy, 0,
                  sizeof(::CollieIterator) );
      instance.SetDelete(&delete_CollieIterator);
      instance.SetDeleteArray(&deleteArray_CollieIterator);
      instance.SetDestructor(&destruct_CollieIterator);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CollieIterator*)
   {
      return GenerateInitInstanceLocal((::CollieIterator*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::CollieIterator*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *CollieIterator_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::CollieIterator*)0x0)->GetClass();
      CollieIterator_TClassManip(theClass);
   return theClass;
   }

   static void CollieIterator_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *CollieQuickKEYS_Dictionary();
   static void CollieQuickKEYS_TClassManip(TClass*);
   static void *new_CollieQuickKEYS(void *p = 0);
   static void *newArray_CollieQuickKEYS(Long_t size, void *p);
   static void delete_CollieQuickKEYS(void *p);
   static void deleteArray_CollieQuickKEYS(void *p);
   static void destruct_CollieQuickKEYS(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CollieQuickKEYS*)
   {
      ::CollieQuickKEYS *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::CollieQuickKEYS));
      static ::ROOT::TGenericClassInfo 
         instance("CollieQuickKEYS", "CollieQuickKEYS.hh", 24,
                  typeid(::CollieQuickKEYS), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &CollieQuickKEYS_Dictionary, isa_proxy, 0,
                  sizeof(::CollieQuickKEYS) );
      instance.SetNew(&new_CollieQuickKEYS);
      instance.SetNewArray(&newArray_CollieQuickKEYS);
      instance.SetDelete(&delete_CollieQuickKEYS);
      instance.SetDeleteArray(&deleteArray_CollieQuickKEYS);
      instance.SetDestructor(&destruct_CollieQuickKEYS);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CollieQuickKEYS*)
   {
      return GenerateInitInstanceLocal((::CollieQuickKEYS*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::CollieQuickKEYS*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *CollieQuickKEYS_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::CollieQuickKEYS*)0x0)->GetClass();
      CollieQuickKEYS_TClassManip(theClass);
   return theClass;
   }

   static void CollieQuickKEYS_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *CollieXsec_Dictionary();
   static void CollieXsec_TClassManip(TClass*);
   static void *new_CollieXsec(void *p = 0);
   static void *newArray_CollieXsec(Long_t size, void *p);
   static void delete_CollieXsec(void *p);
   static void deleteArray_CollieXsec(void *p);
   static void destruct_CollieXsec(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CollieXsec*)
   {
      ::CollieXsec *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::CollieXsec));
      static ::ROOT::TGenericClassInfo 
         instance("CollieXsec", "CollieXsec.hh", 13,
                  typeid(::CollieXsec), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &CollieXsec_Dictionary, isa_proxy, 0,
                  sizeof(::CollieXsec) );
      instance.SetNew(&new_CollieXsec);
      instance.SetNewArray(&newArray_CollieXsec);
      instance.SetDelete(&delete_CollieXsec);
      instance.SetDeleteArray(&deleteArray_CollieXsec);
      instance.SetDestructor(&destruct_CollieXsec);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CollieXsec*)
   {
      return GenerateInitInstanceLocal((::CollieXsec*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::CollieXsec*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *CollieXsec_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::CollieXsec*)0x0)->GetClass();
      CollieXsec_TClassManip(theClass);
   return theClass;
   }

   static void CollieXsec_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *CollieIOFile_Dictionary();
   static void CollieIOFile_TClassManip(TClass*);
   static void *new_CollieIOFile(void *p = 0);
   static void *newArray_CollieIOFile(Long_t size, void *p);
   static void delete_CollieIOFile(void *p);
   static void deleteArray_CollieIOFile(void *p);
   static void destruct_CollieIOFile(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CollieIOFile*)
   {
      ::CollieIOFile *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::CollieIOFile));
      static ::ROOT::TGenericClassInfo 
         instance("CollieIOFile", "CollieIOFile.hh", 15,
                  typeid(::CollieIOFile), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &CollieIOFile_Dictionary, isa_proxy, 0,
                  sizeof(::CollieIOFile) );
      instance.SetNew(&new_CollieIOFile);
      instance.SetNewArray(&newArray_CollieIOFile);
      instance.SetDelete(&delete_CollieIOFile);
      instance.SetDeleteArray(&deleteArray_CollieIOFile);
      instance.SetDestructor(&destruct_CollieIOFile);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CollieIOFile*)
   {
      return GenerateInitInstanceLocal((::CollieIOFile*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::CollieIOFile*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *CollieIOFile_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::CollieIOFile*)0x0)->GetClass();
      CollieIOFile_TClassManip(theClass);
   return theClass;
   }

   static void CollieIOFile_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_CDF_DZero_Distribution(void *p = 0);
   static void *newArray_CDF_DZero_Distribution(Long_t size, void *p);
   static void delete_CDF_DZero_Distribution(void *p);
   static void deleteArray_CDF_DZero_Distribution(void *p);
   static void destruct_CDF_DZero_Distribution(void *p);
   static void streamer_CDF_DZero_Distribution(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CDF_DZero_Distribution*)
   {
      ::CDF_DZero_Distribution *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::CDF_DZero_Distribution >(0);
      static ::ROOT::TGenericClassInfo 
         instance("CDF_DZero_Distribution", ::CDF_DZero_Distribution::Class_Version(), "CDF_DZero_Distribution.hh", 14,
                  typeid(::CDF_DZero_Distribution), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::CDF_DZero_Distribution::Dictionary, isa_proxy, 16,
                  sizeof(::CDF_DZero_Distribution) );
      instance.SetNew(&new_CDF_DZero_Distribution);
      instance.SetNewArray(&newArray_CDF_DZero_Distribution);
      instance.SetDelete(&delete_CDF_DZero_Distribution);
      instance.SetDeleteArray(&deleteArray_CDF_DZero_Distribution);
      instance.SetDestructor(&destruct_CDF_DZero_Distribution);
      instance.SetStreamerFunc(&streamer_CDF_DZero_Distribution);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CDF_DZero_Distribution*)
   {
      return GenerateInitInstanceLocal((::CDF_DZero_Distribution*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::CDF_DZero_Distribution*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_CDF_DZero_IOpoint(void *p = 0);
   static void *newArray_CDF_DZero_IOpoint(Long_t size, void *p);
   static void delete_CDF_DZero_IOpoint(void *p);
   static void deleteArray_CDF_DZero_IOpoint(void *p);
   static void destruct_CDF_DZero_IOpoint(void *p);
   static void streamer_CDF_DZero_IOpoint(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CDF_DZero_IOpoint*)
   {
      ::CDF_DZero_IOpoint *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::CDF_DZero_IOpoint >(0);
      static ::ROOT::TGenericClassInfo 
         instance("CDF_DZero_IOpoint", ::CDF_DZero_IOpoint::Class_Version(), "CDF_DZero_IOpoint.hh", 15,
                  typeid(::CDF_DZero_IOpoint), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::CDF_DZero_IOpoint::Dictionary, isa_proxy, 16,
                  sizeof(::CDF_DZero_IOpoint) );
      instance.SetNew(&new_CDF_DZero_IOpoint);
      instance.SetNewArray(&newArray_CDF_DZero_IOpoint);
      instance.SetDelete(&delete_CDF_DZero_IOpoint);
      instance.SetDeleteArray(&deleteArray_CDF_DZero_IOpoint);
      instance.SetDestructor(&destruct_CDF_DZero_IOpoint);
      instance.SetStreamerFunc(&streamer_CDF_DZero_IOpoint);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CDF_DZero_IOpoint*)
   {
      return GenerateInitInstanceLocal((::CDF_DZero_IOpoint*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::CDF_DZero_IOpoint*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_CDF_DZero_IOfile(void *p = 0);
   static void *newArray_CDF_DZero_IOfile(Long_t size, void *p);
   static void delete_CDF_DZero_IOfile(void *p);
   static void deleteArray_CDF_DZero_IOfile(void *p);
   static void destruct_CDF_DZero_IOfile(void *p);
   static void streamer_CDF_DZero_IOfile(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CDF_DZero_IOfile*)
   {
      ::CDF_DZero_IOfile *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::CDF_DZero_IOfile >(0);
      static ::ROOT::TGenericClassInfo 
         instance("CDF_DZero_IOfile", ::CDF_DZero_IOfile::Class_Version(), "CDF_DZero_IOfile.hh", 11,
                  typeid(::CDF_DZero_IOfile), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::CDF_DZero_IOfile::Dictionary, isa_proxy, 16,
                  sizeof(::CDF_DZero_IOfile) );
      instance.SetNew(&new_CDF_DZero_IOfile);
      instance.SetNewArray(&newArray_CDF_DZero_IOfile);
      instance.SetDelete(&delete_CDF_DZero_IOfile);
      instance.SetDeleteArray(&deleteArray_CDF_DZero_IOfile);
      instance.SetDestructor(&destruct_CDF_DZero_IOfile);
      instance.SetStreamerFunc(&streamer_CDF_DZero_IOfile);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CDF_DZero_IOfile*)
   {
      return GenerateInitInstanceLocal((::CDF_DZero_IOfile*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::CDF_DZero_IOfile*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr InfoForEfficiencyCalculation::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *InfoForEfficiencyCalculation::Class_Name()
{
   return "InfoForEfficiencyCalculation";
}

//______________________________________________________________________________
const char *InfoForEfficiencyCalculation::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::InfoForEfficiencyCalculation*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int InfoForEfficiencyCalculation::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::InfoForEfficiencyCalculation*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *InfoForEfficiencyCalculation::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::InfoForEfficiencyCalculation*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *InfoForEfficiencyCalculation::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::InfoForEfficiencyCalculation*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr CollieDistribution::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *CollieDistribution::Class_Name()
{
   return "CollieDistribution";
}

//______________________________________________________________________________
const char *CollieDistribution::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CollieDistribution*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int CollieDistribution::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CollieDistribution*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *CollieDistribution::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CollieDistribution*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *CollieDistribution::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CollieDistribution*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr CollieEventList::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *CollieEventList::Class_Name()
{
   return "CollieEventList";
}

//______________________________________________________________________________
const char *CollieEventList::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CollieEventList*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int CollieEventList::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CollieEventList*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *CollieEventList::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CollieEventList*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *CollieEventList::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CollieEventList*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr CollieMasspoint::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *CollieMasspoint::Class_Name()
{
   return "CollieMasspoint";
}

//______________________________________________________________________________
const char *CollieMasspoint::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CollieMasspoint*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int CollieMasspoint::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CollieMasspoint*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *CollieMasspoint::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CollieMasspoint*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *CollieMasspoint::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CollieMasspoint*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr CollieChannel::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *CollieChannel::Class_Name()
{
   return "CollieChannel";
}

//______________________________________________________________________________
const char *CollieChannel::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CollieChannel*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int CollieChannel::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CollieChannel*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *CollieChannel::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CollieChannel*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *CollieChannel::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CollieChannel*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr CDF_DZero_Distribution::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *CDF_DZero_Distribution::Class_Name()
{
   return "CDF_DZero_Distribution";
}

//______________________________________________________________________________
const char *CDF_DZero_Distribution::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CDF_DZero_Distribution*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int CDF_DZero_Distribution::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CDF_DZero_Distribution*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *CDF_DZero_Distribution::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CDF_DZero_Distribution*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *CDF_DZero_Distribution::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CDF_DZero_Distribution*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr CDF_DZero_IOpoint::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *CDF_DZero_IOpoint::Class_Name()
{
   return "CDF_DZero_IOpoint";
}

//______________________________________________________________________________
const char *CDF_DZero_IOpoint::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CDF_DZero_IOpoint*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int CDF_DZero_IOpoint::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CDF_DZero_IOpoint*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *CDF_DZero_IOpoint::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CDF_DZero_IOpoint*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *CDF_DZero_IOpoint::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CDF_DZero_IOpoint*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr CDF_DZero_IOfile::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *CDF_DZero_IOfile::Class_Name()
{
   return "CDF_DZero_IOfile";
}

//______________________________________________________________________________
const char *CDF_DZero_IOfile::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CDF_DZero_IOfile*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int CDF_DZero_IOfile::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::CDF_DZero_IOfile*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *CDF_DZero_IOfile::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CDF_DZero_IOfile*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *CDF_DZero_IOfile::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::CDF_DZero_IOfile*)0x0)->GetClass(); }
   return fgIsA;
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_CollieHistogram(void *p) {
      return  p ? new(p) ::CollieHistogram : new ::CollieHistogram;
   }
   static void *newArray_CollieHistogram(Long_t nElements, void *p) {
      return p ? new(p) ::CollieHistogram[nElements] : new ::CollieHistogram[nElements];
   }
   // Wrapper around operator delete
   static void delete_CollieHistogram(void *p) {
      delete ((::CollieHistogram*)p);
   }
   static void deleteArray_CollieHistogram(void *p) {
      delete [] ((::CollieHistogram*)p);
   }
   static void destruct_CollieHistogram(void *p) {
      typedef ::CollieHistogram current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::CollieHistogram

//______________________________________________________________________________
void InfoForEfficiencyCalculation::Streamer(TBuffer &R__b)
{
   // Stream an object of class InfoForEfficiencyCalculation.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      R__b >> sigmaP;
      R__b >> sigmaN;
      R__b >> assym;
      R__b >> exclusionSum;
      R__b >> baseDeltaEfficiency;
      R__b.CheckByteCount(R__s, R__c, InfoForEfficiencyCalculation::IsA());
   } else {
      R__c = R__b.WriteVersion(InfoForEfficiencyCalculation::IsA(), kTRUE);
      R__b << sigmaP;
      R__b << sigmaN;
      R__b << assym;
      R__b << exclusionSum;
      R__b << baseDeltaEfficiency;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_InfoForEfficiencyCalculation(void *p) {
      return  p ? new(p) ::InfoForEfficiencyCalculation : new ::InfoForEfficiencyCalculation;
   }
   static void *newArray_InfoForEfficiencyCalculation(Long_t nElements, void *p) {
      return p ? new(p) ::InfoForEfficiencyCalculation[nElements] : new ::InfoForEfficiencyCalculation[nElements];
   }
   // Wrapper around operator delete
   static void delete_InfoForEfficiencyCalculation(void *p) {
      delete ((::InfoForEfficiencyCalculation*)p);
   }
   static void deleteArray_InfoForEfficiencyCalculation(void *p) {
      delete [] ((::InfoForEfficiencyCalculation*)p);
   }
   static void destruct_InfoForEfficiencyCalculation(void *p) {
      typedef ::InfoForEfficiencyCalculation current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_InfoForEfficiencyCalculation(TBuffer &buf, void *obj) {
      ((::InfoForEfficiencyCalculation*)obj)->::InfoForEfficiencyCalculation::Streamer(buf);
   }
} // end of namespace ROOT for class ::InfoForEfficiencyCalculation

//______________________________________________________________________________
void CollieDistribution::Streamer(TBuffer &R__b)
{
   // Stream an object of class CollieDistribution.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TNamed::Streamer(R__b);
      R__b >> fMinX;
      R__b >> fMaxX;
      R__b >> fMinY;
      R__b >> fMaxY;
      R__b >> fEfficiency;
      R__b >> fNXbins;
      R__b >> fNYbins;
      R__b >> fTrueBinCount;
      fBins.Streamer(R__b);
      fBinStat.Streamer(R__b);
      R__b >> fNmodels;
      fModelXsecs.Streamer(R__b);
      fSystNames.Streamer(R__b);
      fSystematicsPos.Streamer(R__b);
      fSystematicsNeg.Streamer(R__b);
      fFloatFlag.Streamer(R__b);
      fLogNormalFlag.Streamer(R__b);
      R__b >> fPoissonFlag;
      R__b >> fPoissonNorm;
      R__b >> fPoissonErrPos;
      R__b >> fPoissonErrNeg;
      R__b >> fNsyst;
      R__b >> fLinearized;
      R__b.CheckByteCount(R__s, R__c, CollieDistribution::IsA());
   } else {
      R__c = R__b.WriteVersion(CollieDistribution::IsA(), kTRUE);
      TNamed::Streamer(R__b);
      R__b << fMinX;
      R__b << fMaxX;
      R__b << fMinY;
      R__b << fMaxY;
      R__b << fEfficiency;
      R__b << fNXbins;
      R__b << fNYbins;
      R__b << fTrueBinCount;
      fBins.Streamer(R__b);
      fBinStat.Streamer(R__b);
      R__b << fNmodels;
      fModelXsecs.Streamer(R__b);
      fSystNames.Streamer(R__b);
      fSystematicsPos.Streamer(R__b);
      fSystematicsNeg.Streamer(R__b);
      fFloatFlag.Streamer(R__b);
      fLogNormalFlag.Streamer(R__b);
      R__b << fPoissonFlag;
      R__b << fPoissonNorm;
      R__b << fPoissonErrPos;
      R__b << fPoissonErrNeg;
      R__b << fNsyst;
      R__b << fLinearized;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_CollieDistribution(void *p) {
      return  p ? new(p) ::CollieDistribution : new ::CollieDistribution;
   }
   static void *newArray_CollieDistribution(Long_t nElements, void *p) {
      return p ? new(p) ::CollieDistribution[nElements] : new ::CollieDistribution[nElements];
   }
   // Wrapper around operator delete
   static void delete_CollieDistribution(void *p) {
      delete ((::CollieDistribution*)p);
   }
   static void deleteArray_CollieDistribution(void *p) {
      delete [] ((::CollieDistribution*)p);
   }
   static void destruct_CollieDistribution(void *p) {
      typedef ::CollieDistribution current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_CollieDistribution(TBuffer &buf, void *obj) {
      ((::CollieDistribution*)obj)->::CollieDistribution::Streamer(buf);
   }
} // end of namespace ROOT for class ::CollieDistribution

//______________________________________________________________________________
void CollieEventList::Streamer(TBuffer &R__b)
{
   // Stream an object of class CollieEventList.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TNamed::Streamer(R__b);
      fName.Streamer(R__b);
      R__b >> fDimensionality;
      R__b >> fNdoubles;
      R__b >> fNints;
      R__b >> fNevents;
      fDoubles.Streamer(R__b);
      fIntegers.Streamer(R__b);
      fValueNames.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, CollieEventList::IsA());
   } else {
      R__c = R__b.WriteVersion(CollieEventList::IsA(), kTRUE);
      TNamed::Streamer(R__b);
      fName.Streamer(R__b);
      R__b << fDimensionality;
      R__b << fNdoubles;
      R__b << fNints;
      R__b << fNevents;
      fDoubles.Streamer(R__b);
      fIntegers.Streamer(R__b);
      fValueNames.Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_CollieEventList(void *p) {
      return  p ? new(p) ::CollieEventList : new ::CollieEventList;
   }
   static void *newArray_CollieEventList(Long_t nElements, void *p) {
      return p ? new(p) ::CollieEventList[nElements] : new ::CollieEventList[nElements];
   }
   // Wrapper around operator delete
   static void delete_CollieEventList(void *p) {
      delete ((::CollieEventList*)p);
   }
   static void deleteArray_CollieEventList(void *p) {
      delete [] ((::CollieEventList*)p);
   }
   static void destruct_CollieEventList(void *p) {
      typedef ::CollieEventList current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_CollieEventList(TBuffer &buf, void *obj) {
      ((::CollieEventList*)obj)->::CollieEventList::Streamer(buf);
   }
} // end of namespace ROOT for class ::CollieEventList

//______________________________________________________________________________
void CollieMasspoint::Streamer(TBuffer &R__b)
{
   // Stream an object of class CollieMasspoint.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TNamed::Streamer(R__b);
      R__b >> fIndepVar1;
      R__b >> fIndepVar2;
      R__b >> fIndepVar3;
      R__b >> fNXbins;
      R__b >> fNYbins;
      R__b >> fSigSyst;
      R__b >> fBkgSyst;
      R__b >> fMinX;
      R__b >> fMaxX;
      R__b >> fMinY;
      R__b >> fMaxY;
      fSignals.Streamer(R__b);
      fBackgrounds.Streamer(R__b);
      R__b >> fDataDistribution;
      R__b >> fDataEventList;
      R__b.CheckByteCount(R__s, R__c, CollieMasspoint::IsA());
   } else {
      R__c = R__b.WriteVersion(CollieMasspoint::IsA(), kTRUE);
      TNamed::Streamer(R__b);
      R__b << fIndepVar1;
      R__b << fIndepVar2;
      R__b << fIndepVar3;
      R__b << fNXbins;
      R__b << fNYbins;
      R__b << fSigSyst;
      R__b << fBkgSyst;
      R__b << fMinX;
      R__b << fMaxX;
      R__b << fMinY;
      R__b << fMaxY;
      fSignals.Streamer(R__b);
      fBackgrounds.Streamer(R__b);
      R__b << fDataDistribution;
      R__b << fDataEventList;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_CollieMasspoint(void *p) {
      return  p ? new(p) ::CollieMasspoint : new ::CollieMasspoint;
   }
   static void *newArray_CollieMasspoint(Long_t nElements, void *p) {
      return p ? new(p) ::CollieMasspoint[nElements] : new ::CollieMasspoint[nElements];
   }
   // Wrapper around operator delete
   static void delete_CollieMasspoint(void *p) {
      delete ((::CollieMasspoint*)p);
   }
   static void deleteArray_CollieMasspoint(void *p) {
      delete [] ((::CollieMasspoint*)p);
   }
   static void destruct_CollieMasspoint(void *p) {
      typedef ::CollieMasspoint current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_CollieMasspoint(TBuffer &buf, void *obj) {
      ((::CollieMasspoint*)obj)->::CollieMasspoint::Streamer(buf);
   }
} // end of namespace ROOT for class ::CollieMasspoint

//______________________________________________________________________________
void CollieChannel::Streamer(TBuffer &R__b)
{
   // Stream an object of class CollieChannel.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TNamed::Streamer(R__b);
      fChannelName.Streamer(R__b);
      fComments.Streamer(R__b);
      fCollieVersion.Streamer(R__b);
      R__b >> fLuminosity;
      R__b >> fNindepVars;
      int R__i;
      for (R__i = 0; R__i < 3; R__i++)
         fIndepVarName[R__i].Streamer(R__b);
      R__b >> fNsignals;
      fSignalNames.Streamer(R__b);
      R__b >> fNbackgrounds;
      fBkgdNames.Streamer(R__b);
      R__b >> fNmodels;
      fModelNames.Streamer(R__b);
      fCreationComputer.Streamer(R__b);
      fCreationUser.Streamer(R__b);
      fCreationTime.Streamer(R__b);
      fMasspoints_var1.Streamer(R__b);
      fMasspoints_var2.Streamer(R__b);
      fMasspoints_var3.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, CollieChannel::IsA());
   } else {
      R__c = R__b.WriteVersion(CollieChannel::IsA(), kTRUE);
      TNamed::Streamer(R__b);
      fChannelName.Streamer(R__b);
      fComments.Streamer(R__b);
      fCollieVersion.Streamer(R__b);
      R__b << fLuminosity;
      R__b << fNindepVars;
      int R__i;
      for (R__i = 0; R__i < 3; R__i++)
         fIndepVarName[R__i].Streamer(R__b);
      R__b << fNsignals;
      fSignalNames.Streamer(R__b);
      R__b << fNbackgrounds;
      fBkgdNames.Streamer(R__b);
      R__b << fNmodels;
      fModelNames.Streamer(R__b);
      fCreationComputer.Streamer(R__b);
      fCreationUser.Streamer(R__b);
      fCreationTime.Streamer(R__b);
      fMasspoints_var1.Streamer(R__b);
      fMasspoints_var2.Streamer(R__b);
      fMasspoints_var3.Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_CollieChannel(void *p) {
      return  p ? new(p) ::CollieChannel : new ::CollieChannel;
   }
   static void *newArray_CollieChannel(Long_t nElements, void *p) {
      return p ? new(p) ::CollieChannel[nElements] : new ::CollieChannel[nElements];
   }
   // Wrapper around operator delete
   static void delete_CollieChannel(void *p) {
      delete ((::CollieChannel*)p);
   }
   static void deleteArray_CollieChannel(void *p) {
      delete [] ((::CollieChannel*)p);
   }
   static void destruct_CollieChannel(void *p) {
      typedef ::CollieChannel current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_CollieChannel(TBuffer &buf, void *obj) {
      ((::CollieChannel*)obj)->::CollieChannel::Streamer(buf);
   }
} // end of namespace ROOT for class ::CollieChannel

namespace ROOT {
   // Wrapper around operator delete
   static void delete_CollieIterator(void *p) {
      delete ((::CollieIterator*)p);
   }
   static void deleteArray_CollieIterator(void *p) {
      delete [] ((::CollieIterator*)p);
   }
   static void destruct_CollieIterator(void *p) {
      typedef ::CollieIterator current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::CollieIterator

namespace ROOT {
   // Wrappers around operator new
   static void *new_CollieQuickKEYS(void *p) {
      return  p ? new(p) ::CollieQuickKEYS : new ::CollieQuickKEYS;
   }
   static void *newArray_CollieQuickKEYS(Long_t nElements, void *p) {
      return p ? new(p) ::CollieQuickKEYS[nElements] : new ::CollieQuickKEYS[nElements];
   }
   // Wrapper around operator delete
   static void delete_CollieQuickKEYS(void *p) {
      delete ((::CollieQuickKEYS*)p);
   }
   static void deleteArray_CollieQuickKEYS(void *p) {
      delete [] ((::CollieQuickKEYS*)p);
   }
   static void destruct_CollieQuickKEYS(void *p) {
      typedef ::CollieQuickKEYS current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::CollieQuickKEYS

namespace ROOT {
   // Wrappers around operator new
   static void *new_CollieXsec(void *p) {
      return  p ? new(p) ::CollieXsec : new ::CollieXsec;
   }
   static void *newArray_CollieXsec(Long_t nElements, void *p) {
      return p ? new(p) ::CollieXsec[nElements] : new ::CollieXsec[nElements];
   }
   // Wrapper around operator delete
   static void delete_CollieXsec(void *p) {
      delete ((::CollieXsec*)p);
   }
   static void deleteArray_CollieXsec(void *p) {
      delete [] ((::CollieXsec*)p);
   }
   static void destruct_CollieXsec(void *p) {
      typedef ::CollieXsec current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::CollieXsec

namespace ROOT {
   // Wrappers around operator new
   static void *new_CollieIOFile(void *p) {
      return  p ? new(p) ::CollieIOFile : new ::CollieIOFile;
   }
   static void *newArray_CollieIOFile(Long_t nElements, void *p) {
      return p ? new(p) ::CollieIOFile[nElements] : new ::CollieIOFile[nElements];
   }
   // Wrapper around operator delete
   static void delete_CollieIOFile(void *p) {
      delete ((::CollieIOFile*)p);
   }
   static void deleteArray_CollieIOFile(void *p) {
      delete [] ((::CollieIOFile*)p);
   }
   static void destruct_CollieIOFile(void *p) {
      typedef ::CollieIOFile current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::CollieIOFile

//______________________________________________________________________________
void CDF_DZero_Distribution::Streamer(TBuffer &R__b)
{
   // Stream an object of class CDF_DZero_Distribution.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TNamed::Streamer(R__b);
      R__b >> fDistribution;
      fSystematicsPos.Streamer(R__b);
      fSystematicsNeg.Streamer(R__b);
      fSystematicsNames.Streamer(R__b);
      R__b >> f2Ddistribution;
      R__b >> fNXbins;
      R__b >> fNYbins;
      R__b >> fMinX;
      R__b >> fMaxX;
      R__b >> fMinY;
      R__b >> fMaxY;
      R__b.CheckByteCount(R__s, R__c, CDF_DZero_Distribution::IsA());
   } else {
      R__c = R__b.WriteVersion(CDF_DZero_Distribution::IsA(), kTRUE);
      TNamed::Streamer(R__b);
      R__b << fDistribution;
      fSystematicsPos.Streamer(R__b);
      fSystematicsNeg.Streamer(R__b);
      fSystematicsNames.Streamer(R__b);
      R__b << f2Ddistribution;
      R__b << fNXbins;
      R__b << fNYbins;
      R__b << fMinX;
      R__b << fMaxX;
      R__b << fMinY;
      R__b << fMaxY;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_CDF_DZero_Distribution(void *p) {
      return  p ? new(p) ::CDF_DZero_Distribution : new ::CDF_DZero_Distribution;
   }
   static void *newArray_CDF_DZero_Distribution(Long_t nElements, void *p) {
      return p ? new(p) ::CDF_DZero_Distribution[nElements] : new ::CDF_DZero_Distribution[nElements];
   }
   // Wrapper around operator delete
   static void delete_CDF_DZero_Distribution(void *p) {
      delete ((::CDF_DZero_Distribution*)p);
   }
   static void deleteArray_CDF_DZero_Distribution(void *p) {
      delete [] ((::CDF_DZero_Distribution*)p);
   }
   static void destruct_CDF_DZero_Distribution(void *p) {
      typedef ::CDF_DZero_Distribution current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_CDF_DZero_Distribution(TBuffer &buf, void *obj) {
      ((::CDF_DZero_Distribution*)obj)->::CDF_DZero_Distribution::Streamer(buf);
   }
} // end of namespace ROOT for class ::CDF_DZero_Distribution

//______________________________________________________________________________
void CDF_DZero_IOpoint::Streamer(TBuffer &R__b)
{
   // Stream an object of class CDF_DZero_IOpoint.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TNamed::Streamer(R__b);
      fSignals.Streamer(R__b);
      fBackgrounds.Streamer(R__b);
      R__b >> fDataDistribution;
      R__b >> fIndepVar1;
      R__b >> fIndepVar2;
      R__b >> fIndepVar3;
      R__b >> fNXbins;
      R__b >> fNYbins;
      R__b >> fMinX;
      R__b >> fMaxX;
      R__b >> fMinY;
      R__b >> fMaxY;
      R__b.CheckByteCount(R__s, R__c, CDF_DZero_IOpoint::IsA());
   } else {
      R__c = R__b.WriteVersion(CDF_DZero_IOpoint::IsA(), kTRUE);
      TNamed::Streamer(R__b);
      fSignals.Streamer(R__b);
      fBackgrounds.Streamer(R__b);
      R__b << fDataDistribution;
      R__b << fIndepVar1;
      R__b << fIndepVar2;
      R__b << fIndepVar3;
      R__b << fNXbins;
      R__b << fNYbins;
      R__b << fMinX;
      R__b << fMaxX;
      R__b << fMinY;
      R__b << fMaxY;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_CDF_DZero_IOpoint(void *p) {
      return  p ? new(p) ::CDF_DZero_IOpoint : new ::CDF_DZero_IOpoint;
   }
   static void *newArray_CDF_DZero_IOpoint(Long_t nElements, void *p) {
      return p ? new(p) ::CDF_DZero_IOpoint[nElements] : new ::CDF_DZero_IOpoint[nElements];
   }
   // Wrapper around operator delete
   static void delete_CDF_DZero_IOpoint(void *p) {
      delete ((::CDF_DZero_IOpoint*)p);
   }
   static void deleteArray_CDF_DZero_IOpoint(void *p) {
      delete [] ((::CDF_DZero_IOpoint*)p);
   }
   static void destruct_CDF_DZero_IOpoint(void *p) {
      typedef ::CDF_DZero_IOpoint current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_CDF_DZero_IOpoint(TBuffer &buf, void *obj) {
      ((::CDF_DZero_IOpoint*)obj)->::CDF_DZero_IOpoint::Streamer(buf);
   }
} // end of namespace ROOT for class ::CDF_DZero_IOpoint

//______________________________________________________________________________
void CDF_DZero_IOfile::Streamer(TBuffer &R__b)
{
   // Stream an object of class CDF_DZero_IOfile.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TNamed::Streamer(R__b);
      fChannelName.Streamer(R__b);
      fCreationComputer.Streamer(R__b);
      fCreationUser.Streamer(R__b);
      fCreationTime.Streamer(R__b);
      R__b >> fNindepVars;
      int R__i;
      for (R__i = 0; R__i < 3; R__i++)
         fIndepVarName[R__i].Streamer(R__b);
      fPoints_var1.Streamer(R__b);
      fPoints_var2.Streamer(R__b);
      fPoints_var3.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, CDF_DZero_IOfile::IsA());
   } else {
      R__c = R__b.WriteVersion(CDF_DZero_IOfile::IsA(), kTRUE);
      TNamed::Streamer(R__b);
      fChannelName.Streamer(R__b);
      fCreationComputer.Streamer(R__b);
      fCreationUser.Streamer(R__b);
      fCreationTime.Streamer(R__b);
      R__b << fNindepVars;
      int R__i;
      for (R__i = 0; R__i < 3; R__i++)
         fIndepVarName[R__i].Streamer(R__b);
      fPoints_var1.Streamer(R__b);
      fPoints_var2.Streamer(R__b);
      fPoints_var3.Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_CDF_DZero_IOfile(void *p) {
      return  p ? new(p) ::CDF_DZero_IOfile : new ::CDF_DZero_IOfile;
   }
   static void *newArray_CDF_DZero_IOfile(Long_t nElements, void *p) {
      return p ? new(p) ::CDF_DZero_IOfile[nElements] : new ::CDF_DZero_IOfile[nElements];
   }
   // Wrapper around operator delete
   static void delete_CDF_DZero_IOfile(void *p) {
      delete ((::CDF_DZero_IOfile*)p);
   }
   static void deleteArray_CDF_DZero_IOfile(void *p) {
      delete [] ((::CDF_DZero_IOfile*)p);
   }
   static void destruct_CDF_DZero_IOfile(void *p) {
      typedef ::CDF_DZero_IOfile current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_CDF_DZero_IOfile(TBuffer &buf, void *obj) {
      ((::CDF_DZero_IOfile*)obj)->::CDF_DZero_IOfile::Streamer(buf);
   }
} // end of namespace ROOT for class ::CDF_DZero_IOfile

namespace {
  void TriggerDictionaryInitialization_collieio_dict_Impl() {
    static const char* headers[] = {
"io/include//CollieDistribution.hh",
"io/include//CollieEventList.hh",
"io/include//CollieChannel.hh",
"io/include//CollieIterator.hh",
"io/include//CollieHistogram.hh",
"io/include//CollieQuickKEYS.hh",
"io/include//CollieIOFile.hh",
"io/include//CollieXsec.hh",
"io/include//CollieMasspoint.hh",
"io/include//CDF_DZero_Distribution.hh",
"io/include//CDF_DZero_IOpoint.hh",
"io/include//CDF_DZero_IOfile.hh",
"io/include//InfoForEfficiencyCalculation.hh",
0
    };
    static const char* includePaths[] = {
"/pc2014-data2/pguzowski/snsw_3.0/CadfaelBrew/Cellar/root6/6.08.06/include/root",
"io/include/",
"limit/include/",
"CLHEP/Random/",
"minuit/include/",
"/pc2014-data2/pguzowski/snsw_3.0/CadfaelBrew/Cellar/root6/6.08.06/include/root",
"/pc2014-data5/nlane/hnl/collie/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "collieio_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$CollieHistogram.hh")))  __attribute__((annotate("$clingAutoload$io/include//CollieDistribution.hh")))  CollieHistogram;
struct __attribute__((annotate("$clingAutoload$InfoForEfficiencyCalculation.hh")))  __attribute__((annotate("$clingAutoload$io/include//CollieDistribution.hh")))  InfoForEfficiencyCalculation;
class __attribute__((annotate("$clingAutoload$io/include//CollieDistribution.hh")))  CollieDistribution;
class __attribute__((annotate("$clingAutoload$io/include//CollieEventList.hh")))  CollieEventList;
class __attribute__((annotate("$clingAutoload$CollieMasspoint.hh")))  __attribute__((annotate("$clingAutoload$io/include//CollieChannel.hh")))  CollieMasspoint;
class __attribute__((annotate("$clingAutoload$io/include//CollieChannel.hh")))  CollieChannel;
class __attribute__((annotate("$clingAutoload$io/include//CollieIterator.hh")))  CollieIterator;
class __attribute__((annotate("$clingAutoload$io/include//CollieQuickKEYS.hh")))  CollieQuickKEYS;
class __attribute__((annotate("$clingAutoload$CollieXsec.hh")))  __attribute__((annotate("$clingAutoload$io/include//CollieIOFile.hh")))  CollieXsec;
class __attribute__((annotate("$clingAutoload$io/include//CollieIOFile.hh")))  CollieIOFile;
class __attribute__((annotate("$clingAutoload$io/include//CDF_DZero_Distribution.hh")))  CDF_DZero_Distribution;
class __attribute__((annotate("$clingAutoload$io/include//CDF_DZero_IOpoint.hh")))  CDF_DZero_IOpoint;
class __attribute__((annotate("$clingAutoload$io/include//CDF_DZero_IOfile.hh")))  CDF_DZero_IOfile;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "collieio_dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "io/include//CollieDistribution.hh"
#include "io/include//CollieEventList.hh"
#include "io/include//CollieChannel.hh"
#include "io/include//CollieIterator.hh"
#include "io/include//CollieHistogram.hh"
#include "io/include//CollieQuickKEYS.hh"
#include "io/include//CollieIOFile.hh"
#include "io/include//CollieXsec.hh"
#include "io/include//CollieMasspoint.hh"
#include "io/include//CDF_DZero_Distribution.hh"
#include "io/include//CDF_DZero_IOpoint.hh"
#include "io/include//CDF_DZero_IOfile.hh"
#include "io/include//InfoForEfficiencyCalculation.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"CDF_DZero_Distribution", payloadCode, "@",
"CDF_DZero_IOfile", payloadCode, "@",
"CDF_DZero_IOpoint", payloadCode, "@",
"CollieChannel", payloadCode, "@",
"CollieDistribution", payloadCode, "@",
"CollieEventList", payloadCode, "@",
"CollieHistogram", payloadCode, "@",
"CollieIOFile", payloadCode, "@",
"CollieIterator", payloadCode, "@",
"CollieMasspoint", payloadCode, "@",
"CollieQuickKEYS", payloadCode, "@",
"CollieXsec", payloadCode, "@",
"InfoForEfficiencyCalculation", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("collieio_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_collieio_dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_collieio_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_collieio_dict() {
  TriggerDictionaryInitialization_collieio_dict_Impl();
}
