# Copyright (c) 2013-2022, The MercuryDPM Developers Team. All rights reserved.
# For the list of developers, see <http://www.MercuryDPM.org/Team>.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#   * Neither the name MercuryDPM nor the
#     names of its contributors may be used to endorse or promote products
#     derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE

# This class is used for colored terminal output.

class colorClass:

    def __init__(self):
        self._use_color = True
        self.def_colors()

    def set_use_colors(self, colors):
        self._use_color = colors

    def def_colors(self):
        if self._use_color:
            self.END      = '\33[0m'
            self.BOLD     = '\33[1m'
            self.ITALIC   = '\33[3m'
            self.URL      = '\33[4m'
            self.BLINK    = '\33[5m'
            self.BLINK2   = '\33[6m'
            self.SELECTED = '\33[7m'

            self.BLACK  = '\33[30m'
            self.RED    = '\33[31m'
            self.GREEN  = '\33[32m'
            self.YELLOW = '\33[33m'
            self.BLUE   = '\33[34m'
            self.VIOLET = '\33[35m'
            self.BEIGE  = '\33[36m'
            self.WHITE  = '\33[37m'

            self.BLACKBG  = '\33[40m'
            self.REDBG    = '\33[41m'
            self.GREENBG  = '\33[42m'
            self.YELLOWBG = '\33[43m'
            self.BLUEBG   = '\33[44m'
            self.VIOLETBG = '\33[45m'
            self.BEIGEBG  = '\33[46m'
            self.WHITEBG  = '\33[47m'

            self.GREY    = '\33[90m'
            self.RED2    = '\33[91m'
            self.GREEN2  = '\33[92m'
            self.YELLOW2 = '\33[93m'
            self.BLUE2   = '\33[94m'
            self.VIOLET2 = '\33[95m'
            self.BEIGE2  = '\33[96m'
            self.WHITE2  = '\33[97m'
        else:
            self.END = ''
            self.BOLD = ''
            self.ITALIC = ''
            self.URL = ''
            self.BLINK = ''
            self.BLINK2 = ''
            self.SELECTED = ''

            self.BLACK = ''
            self.RED = ''
            self.GREEN = ''
            self.YELLOW = ''
            self.BLUE = ''
            self.VIOLET = ''
            self.BEIGE = ''
            self.WHITE = ''

            self.BLACKBG = ''
            self.REDBG = ''
            self.GREENBG = ''
            self.YELLOWBG = ''
            self.BLUEBG = ''
            self.VIOLETBG = ''
            self.BEIGEBG = ''
            self.WHITEBG = ''

            self.GREY = ''
            self.RED2 = ''
            self.GREEN2 = ''
            self.YELLOW2 = ''
            self.BLUE2 = ''
            self.VIOLET2 = ''
            self.BEIGE2 = ''
            self.WHITE2 = ''
