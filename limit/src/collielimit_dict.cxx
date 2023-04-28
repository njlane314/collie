// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME limitdIsrcdIdIcollielimit_dict

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
#include "limit/include//CLpoint.hh"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *CLpoint_Dictionary();
   static void CLpoint_TClassManip(TClass*);
   static void *new_CLpoint(void *p = 0);
   static void *newArray_CLpoint(Long_t size, void *p);
   static void delete_CLpoint(void *p);
   static void deleteArray_CLpoint(void *p);
   static void destruct_CLpoint(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CLpoint*)
   {
      ::CLpoint *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::CLpoint));
      static ::ROOT::TGenericClassInfo 
         instance("CLpoint", "CLpoint.hh", 10,
                  typeid(::CLpoint), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &CLpoint_Dictionary, isa_proxy, 0,
                  sizeof(::CLpoint) );
      instance.SetNew(&new_CLpoint);
      instance.SetNewArray(&newArray_CLpoint);
      instance.SetDelete(&delete_CLpoint);
      instance.SetDeleteArray(&deleteArray_CLpoint);
      instance.SetDestructor(&destruct_CLpoint);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CLpoint*)
   {
      return GenerateInitInstanceLocal((::CLpoint*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::CLpoint*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *CLpoint_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::CLpoint*)0x0)->GetClass();
      CLpoint_TClassManip(theClass);
   return theClass;
   }

   static void CLpoint_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_CLpoint(void *p) {
      return  p ? new(p) ::CLpoint : new ::CLpoint;
   }
   static void *newArray_CLpoint(Long_t nElements, void *p) {
      return p ? new(p) ::CLpoint[nElements] : new ::CLpoint[nElements];
   }
   // Wrapper around operator delete
   static void delete_CLpoint(void *p) {
      delete ((::CLpoint*)p);
   }
   static void deleteArray_CLpoint(void *p) {
      delete [] ((::CLpoint*)p);
   }
   static void destruct_CLpoint(void *p) {
      typedef ::CLpoint current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::CLpoint

namespace {
  void TriggerDictionaryInitialization_collielimit_dict_Impl() {
    static const char* headers[] = {
"limit/include//CLpoint.hh",
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
#line 1 "collielimit_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
struct __attribute__((annotate("$clingAutoload$limit/include//CLpoint.hh")))  CLpoint;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "collielimit_dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "limit/include//CLpoint.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"CLpoint", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("collielimit_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_collielimit_dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_collielimit_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_collielimit_dict() {
  TriggerDictionaryInitialization_collielimit_dict_Impl();
}
