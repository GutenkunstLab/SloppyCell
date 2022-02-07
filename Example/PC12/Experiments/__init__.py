from __future__ import absolute_import
import SloppyCell.Collections

from . import ErkMekTraverse2EGF as ErkMekTraverse2EGF
from . import Raf1LandrethEGF as Raf1LandrethEGF
from . import Erk1HYaoEGF as Erk1HYaoEGF
from . import ErkMekTraverse2NGF as ErkMekTraverse2NGF
from . import RasGreen1NGF as RasGreen1NGF
from . import Rap1YorkNGF as Rap1YorkNGF
from . import Erk1HYaoNGF as Erk1HYaoNGF
from . import RasGreen1EGF as RasGreen1EGF
from . import Raf1BRafLandrethNGF as Raf1BRafLandrethNGF
from . import ErkMekTraverseHERover as ErkMekTraverseHERover

exptColl = SloppyCell.Collections.\
        ExperimentCollection([
                              ErkMekTraverse2EGF.expt,
                              Raf1LandrethEGF.expt,
                              Erk1HYaoEGF.expt,
                              ErkMekTraverse2NGF.expt,
                              RasGreen1NGF.expt,
                              Rap1YorkNGF.expt,
                              Erk1HYaoNGF.expt,
                              RasGreen1EGF.expt,
                              Raf1BRafLandrethNGF.expt,
                              ErkMekTraverseHERover.expt,
                              ])
