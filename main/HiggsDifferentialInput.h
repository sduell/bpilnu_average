/*
 *  HiggsCombiner: Florian Bernlochner
 */

#ifndef _HiggsInput
#define _HiggsInput

#include "../utils/Utils.h"
#include "HiggsDifferentialCommon.h"

class HiggsInput {

   public:

      // Constructors
      HiggsInput();
      HiggsInput(settings set);
      
      vector <input> GetInput() { return _v_input; };
    
      TEnv *GetTEnv() { return var_input; };
      
            
   private:
   
      void FillInputVector();
   
      settings _set; 
      
      TEnv *var_input;
      
      vector <input> _v_input;


};

#endif
