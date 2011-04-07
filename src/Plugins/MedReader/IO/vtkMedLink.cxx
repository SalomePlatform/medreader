#include "vtkMedLink.h"

#include "vtkObjectFactory.h"

#include "vtkMedUtilities.h"
#include "vtkMedString.h"

#include <string>
using namespace std;

vtkCxxRevisionMacro(vtkMedLink, "$Revision$")
vtkStandardNewMacro(vtkMedLink)

vtkMedLink::vtkMedLink()
{
	this->MedIterator = -1;
	this->MeshName = vtkMedString::New();
	this->MeshName->SetSize(MED_NAME_SIZE);
	this->Link = NULL;
}

vtkMedLink::~vtkMedLink()
{
	this->MeshName->Delete();
	this->SetLink(NULL);
}

const char* vtkMedLink::GetFullLink(const char* originalFileName)
{
#ifdef _WIN32
  static const char sep = '\\';
#else
  static const char sep = '/';
#endif

  if(this->Link == NULL)
    {
    return NULL;
    }

  // First test if the Link is a full path, then return it.
  if(this->Link != NULL && this->Link[0] == sep)
    {
    return this->Link;
    }

  string name = string(originalFileName);
  size_t pos = name.find_last_of(sep);
  if(pos == string::npos)
    {
    return this->Link;
    }

  string clean_link = this->Link;
  string to_remove = string(".") + sep;
  int to_remove_size = to_remove.size();
  while(clean_link.substr(0, to_remove_size) == to_remove)
    clean_link = clean_link.substr(to_remove_size, string::npos);

  string path = name.substr(0, pos+1);
  this->FullLinkPath = path + clean_link;
  return this->FullLinkPath.c_str();
}

void   vtkMedLink::SetMountedIterator(med_class what, med_int mit)
{
  this->MountedIterator[what] = mit;
}

med_int  vtkMedLink::GetMountedIterator(med_class what)
{
  if(this->MountedIterator.find(what) == this->MountedIterator.end())
    return -1;

  return this->MountedIterator[what];
}

void vtkMedLink::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	PRINT_IVAR(os, indent, MedIterator);
}
