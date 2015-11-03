#include "Variant.h"
#include "../fastahack/src/Fasta.h"

int main()
{
  std::string vcfFileName = "../tests/example.vcf";
  std::string fastaFileName = "../tests/reference.fa";

  vcflib::VariantCallFile variantFile;
  variantFile.open(vcfFileName);

  if (!variantFile.is_open())
  {
    std::cerr << "Could not open VCF file." << std::endl;
    exit(1);
  }

  FastaReference reference;

  if (fastaFileName.empty())
  {
    std::cerr << "Invalid reference." << std::endl;
    exit(1);
  }
  else
  {
    reference.open(fastaFileName);
  }

  std::cout << variantFile.header << std::endl;
  long int last_end = 1;
  std::string sequenceName;
  vcflib::Variant var(variantFile);

  while (variantFile.getNextVariant(var))
  {
    if (sequenceName.empty())
    {
      sequenceName = var.sequenceName;
    }
    else if (sequenceName != var.sequenceName)
    {
      // emit last record from previous chrom
      // these should be refactored.....
      vcflib::Variant refvar(variantFile);
      if (var.position - last_end > 0)
      {
        refvar.ref = reference.getSubSequence(sequenceName, last_end - 1, var.position - last_end);
        refvar.quality = 0;
        refvar.position = last_end;
        refvar.sequenceName = sequenceName;
        std::cout << refvar << std::endl;
      }

      last_end = 1;
      sequenceName = var.sequenceName;
    }

    // generate the last reference record if we have sequence between variants
    if (var.position - last_end > 0)
    {
      vcflib::Variant refvar(variantFile);
      refvar.quality = 0;
      refvar.position = last_end;
      refvar.sequenceName = sequenceName;
      refvar.ref = reference.getSubSequence(sequenceName, last_end - 1, var.position - last_end);
      std::cout << refvar << std::endl;
    }

    std::cout << var << std::endl;
    last_end = var.position + var.ref.size();
  }

  if (reference.sequenceLength(sequenceName) - last_end > 0)
  {
    vcflib::Variant refvar(variantFile);
    refvar.quality = 0;
    refvar.position = last_end;
    refvar.sequenceName = sequenceName;
    refvar.ref = reference.getSubSequence(sequenceName, last_end, reference.sequenceLength(sequenceName) - last_end);
    std::cout << refvar << std::endl;
  }
}
