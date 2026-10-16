"""MolFromSequence / MolFromFASTA / MolToSequence / MolToFASTA. Also valid on upstream RDKit."""

import pytest
from rdoxide import Chem


def test_peptide_roundtrip():
    m = Chem.MolFromSequence("AGKV")
    assert Chem.MolToSequence(m) == "AGKV"
    from rdoxide.Chem import rdMolDescriptors
    assert rdMolDescriptors.CalcMolFormula(m) == "C16H31N5O5"
    info = m.GetAtomWithIdx(1).GetPDBResidueInfo()
    assert (info.GetName(), info.GetResidueName(), info.GetChainId()) == (" CA ", "ALA", "A")


def test_d_amino_acids():
    m = Chem.MolFromSequence("aG", flavor=1)
    assert Chem.MolToSequence(m) == "aG"
    assert Chem.MolFromSequence("aG").GetAtomWithIdx(1).GetPDBResidueInfo().GetResidueName() == "ALA"


def test_invalid_characters():
    assert Chem.MolFromSequence("AZ") is None
    assert Chem.MolFromSequence("A", flavor=42) is None


def test_chain_breaks_and_whitespace():
    m = Chem.MolFromSequence("GG.AA")
    chains = {a.GetPDBResidueInfo().GetChainId() for a in m.GetAtoms()}
    assert chains == {"A", "B"}
    # a space ends the chain without starting a new one
    assert len(Chem.GetMolFrags(Chem.MolFromSequence("G G"))) == 2


def test_nucleic_acids():
    rna = Chem.MolFromSequence("ACGU", flavor=2)
    assert Chem.MolToSequence(rna) == "ACGU"
    dna = Chem.MolFromSequence("ACGT", flavor=6)
    assert Chem.MolToSequence(dna) == "ACGT"
    capped = Chem.MolFromSequence("AC", flavor=5)
    plain = Chem.MolFromSequence("AC", flavor=2)
    assert capped.GetNumAtoms() > plain.GetNumAtoms()


def test_fasta():
    m = Chem.MolFromFASTA(">my protein\nACD\n")
    assert m.GetProp("_Name") == "my protein"
    assert Chem.MolToFASTA(m) == ">my protein\nACD\n"
    assert Chem.MolToFASTA(Chem.MolFromSmiles("CC")) == ""


def test_sanitize_flag():
    m = Chem.MolFromSequence("A", sanitize=False)
    assert m.GetNumAtoms() == 6


def test_helm_peptide_roundtrip():
    for helm in ("PEPTIDE1{A.C.D}$$$$", "PEPTIDE1{A.[dA].G}$$$$", "PEPTIDE1{[ac].A.G.[am]}$$$$"):
        m = Chem.MolFromHELM(helm)
        assert m is not None, helm
        assert Chem.MolToHELM(m) == helm


def test_helm_connections():
    helm = "PEPTIDE1{A.C.D}|PEPTIDE2{G.C}$PEPTIDE1,PEPTIDE2,2:R3-2:R3$$$"
    m = Chem.MolFromHELM(helm)
    assert Chem.MolToHELM(m) == helm
    assert len(Chem.GetMolFrags(m)) == 1
    # a head-to-tail cycle replaces the terminal OXT
    cyc = Chem.MolFromHELM("PEPTIDE1{A.C.D.G}$PEPTIDE1,PEPTIDE1,1:R1-4:R2$$$")
    assert cyc.GetNumAtoms() == Chem.MolFromHELM("PEPTIDE1{A.C.D.G}$$$$").GetNumAtoms() - 1


def test_helm_nucleic_acids():
    for helm in ("RNA1{R(A)P.R(C)P.R(G)}$$$$", "RNA1{P.R(A)P.R(C)P}$$$$"):
        m = Chem.MolFromHELM(helm)
        assert Chem.MolToHELM(m) == helm


def test_helm_invalid():
    for bad in ("", "PEPTIDE1{A.X}$$$$", "PEPTIDE1{A.C}$PEPTIDE1,PEPTIDE1,1:R3-2:R3$$$"):
        assert Chem.MolFromHELM(bad) is None
    assert Chem.MolFromHELM("$$$$").GetNumAtoms() == 0
    assert Chem.MolToHELM(Chem.MolFromSmiles("CC")) == ""
