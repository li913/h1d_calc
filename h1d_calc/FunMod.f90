!DIR$ define Van
module FunMod
!DIR$ if defined (Van)
use FunMod_van
!DIR$ else
use FunMod_Kosugi, only: ModelNum, FS, FH, FK, FQ, FC, FKcap
!DIR$ end if 
end module FunMod