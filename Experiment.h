#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <string>
#include "XML.h"

class Experiment {
 public:
  Experiment(const char *filename);
  ~Experiment() {};

  int numSteps() const;
  bool populateNode(XMLNode *root, int trial_idx) const;
 private:

  std::string getIntegerString(double alpha) const;
  std::string getDoubleString(double alpha) const;
  std::string getVectorString(double alpha) const;
  std::string getFileString(double alpha) const;
    
  XMLNode m_experiment;
};

#endif  
