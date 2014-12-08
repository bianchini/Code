#ifndef StandardIncludes_h
#define StandardIncludes_h

#include "interface/Utils.h"

using namespace std;

typedef TLorentzVector LV;

namespace Algo {
   
  class TransferFunction  {    
  public:
    TransferFunction(const string&, const string&, const int);
    ~TransferFunction();
    void         init(const double*);
    const string getFormula() const; 
    double       eval       (const double& , const double&) const ;
    void         add_pdf_obs(const string&, const Object&, const int&);    
    double       get_threshold() const;
    void         set_threshold(const double&) ;
  private:
    map<string,double> pdfs;
    double    get_pdfs() const;
    string    formula;
    TFormula* f;
    double    threshold;
    int       verbose;
  };
  
  class DecayBuilder {
  public:
    virtual ~DecayBuilder() {};
    virtual double eval( const double* , LV&) = 0; 
    virtual void   print(ostream&) = 0;
    virtual Decay  get_decay() = 0;
  };
     
  class CombBuilder: public DecayBuilder {
  public:
    CombBuilder(vector<DecayBuilder*>&& = vector<DecayBuilder*>(), int =0);
    ~CombBuilder();
    void          add(DecayBuilder*);
    double        eval ( const double* , LV&);
    double        eval ( const double* );
    void          print(ostream&);
    DecayBuilder* at(const size_t&);
    size_t        size();
    Decay         get_decay();
  private:
    vector<DecayBuilder*> combined;
    int verbose;
  };

  class METBuilder: public DecayBuilder {
  public:
    METBuilder(int =0);
    ~METBuilder();
    void   init(const Object&);
    double eval (  const double* , LV& );
    void   fix_vars();
    void   print(ostream&);   
    Decay  get_decay();
  private:
    LV p4_invisible;
    TransferFunction* tf_met;
    Decay  decay;
    int    verbose;
    int    saturate;
  };
  
  
  class TopHadBuilder: public DecayBuilder {
    
  public:
    TopHadBuilder(int =0);
    ~TopHadBuilder();
    void   init(const FinalState& , const Object&, const size_t&);
    double eval ( const double* , LV&) ;
    void   print(ostream&);     
    Decay  get_decay();
    vector<size_t> get_variables();
  private:    
    LV p4_q;
    LV p4_qbar;
    LV p4_b;
    size_t index_q;
    size_t index_qbar;
    size_t index_b;
    TransferFunction* tf_q;
    TransferFunction* tf_qbar;
    TransferFunction* tf_b;
    Decay  decay;
    size_t errFlag;
    int    verbose;
    int    qbar_lost;
  };

  class WHadBuilder: public DecayBuilder {    
  public:
    WHadBuilder(int =0);
    ~WHadBuilder();
    void   init(const FinalState& , const Object&, const size_t&);
    double eval ( const double* , LV&) ;
    void   print(ostream&);     
    Decay  get_decay();
    vector<size_t> get_variables();
  private:
    LV p4_q;
    LV p4_qbar;
    size_t index_q;
    size_t index_qbar;
    TransferFunction* tf_q;
    TransferFunction* tf_qbar;
    Decay  decay;
    size_t errFlag;
    int    verbose;
  };

 class TopLepBuilder: public DecayBuilder {    
  public:
    TopLepBuilder(int =0);
    ~TopLepBuilder();
    void   init(const FinalState& , const Object&, const size_t&);
    double eval ( const double* , LV&) ;
    void   print(ostream&);     
    Decay  get_decay();
    vector<size_t> get_variables();
  private:    
    LV p4_l;
    LV p4_b;
    size_t index_l;
    size_t index_b;
    TransferFunction* tf_b;
    Decay  decay;
    size_t errFlag;
    int    verbose;
  };
  
  class RadiationBuilder : public DecayBuilder {    
  public:
    RadiationBuilder(int =0, Algo::Decay =Decay::Radiation_g);
    ~RadiationBuilder();
    void   init(const FinalState& , const Object&, const size_t&);
    double eval ( const double* , LV&) ;    
    void   print(ostream&);     
    Decay  get_decay();
    vector<size_t> get_variables();
  private:
    LV p4_g;
    size_t index_g;
    TransferFunction* tf_g;
    Decay decay;
    int   verbose;
  };
 
  class HiggsBuilder: public DecayBuilder {    
  public:
    HiggsBuilder(int =0);
    ~HiggsBuilder();
    void   init(const FinalState& , const Object&, const size_t&);
    double eval ( const double* , LV&) ;
    void   print(ostream&);     
    Decay  get_decay();
    vector<size_t> get_variables();
  private:
    LV p4_b;
    LV p4_bbar;
    size_t index_b;
    size_t index_bbar;
    TransferFunction* tf_b;
    TransferFunction* tf_bbar;
    Decay  decay;
    size_t errFlag;
    int    verbose;
  };

}

#endif
