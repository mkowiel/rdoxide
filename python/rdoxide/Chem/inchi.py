#
#  Copyright (c) 2011, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# rdoxide: PROVISIONAL version of upstream External/INCHI-API/python/inchi.py
# (which RDKit installs as rdkit/Chem/inchi.py). The InChI port only covers
# InChIKeys so far (docs/INCHI_PORT.md), so this module behaves like RDKit's
# build without InChI support (`noinchi.py`: INCHI_AVAILABLE = False, no
# MolToInchi/MolFromInchi, which the upstream tests use to skip) and adds the
# two functions that already work. It is replaced by the verbatim upstream
# file once MolToInchi and MolFromInchi are ported. The native
# rdoxide.Chem.rdinchi module already has the full upstream surface; its
# unported functions raise NotImplementedError.

INCHI_AVAILABLE = False

from rdoxide.Chem import rdinchi


class InchiReadWriteError(Exception):
  pass


def InchiToInchiKey(inchi):
  """Return the InChI key for the given InChI string. Return None on error"""
  ret = rdinchi.InchiToInchiKey(inchi)
  if ret:
    return ret
  else:
    return None


GetInchiVersion = rdinchi.GetInchiVersion

__all__ = ['InchiReadWriteError', 'InchiToInchiKey', 'GetInchiVersion', 'INCHI_AVAILABLE']
