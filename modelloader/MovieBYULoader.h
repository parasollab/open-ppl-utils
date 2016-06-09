#ifndef MOVIE_BYU_LOADER_H_
#define MOVIE_BYU_LOADER_H_

#include <fstream>
#include <utility>
#include <vector>

#include "IModel.h"

////////////////////////////////////////////////////////////////////////////////
/// \brief A loader for Movie.BYU files.
////////////////////////////////////////////////////////////////////////////////
class CMovieBYULoader : public IModel {

  public:

    ///\name IModel Overriedes
    ///@{

    virtual bool ParseFile(bool _silent = false) override;

    ///@}

  private:

    bool ParseSection1(std::ifstream& _in);
    bool ParseSection2(std::ifstream& _in);
    bool ParseSection3(std::ifstream& _in);

    int m_partsSize{0};
    int m_vertexSize{0};
    int m_polygonSize{0};
    int m_edgeSize{0};

    std::vector<std::pair<int,int>> m_parts;
};

#endif
