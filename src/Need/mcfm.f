      program mcfm
      implicit none
      character*72 inputfile,workdir
      call determinefilenames(inputfile,workdir)
      call mcfmmain(inputfile,workdir)
      end
