package main

import (
	"errors"
	"glycodist/pdb"
	"glycodist/uniprot"
	"math"
	"sort"
)

type residueDistance struct {
	UnpPos   int64
	Residue  *pdb.Residue
	Distance float64
}

func distance(atom1 *pdb.Atom, atom2 *pdb.Atom) float64 {
	return math.Sqrt(math.Pow(atom1.X-atom2.X, 2) + math.Pow(atom1.Y-atom2.Y, 2) + math.Pow(atom1.Z-atom2.Z, 2))
}

func residuesDistance(res1 *pdb.Residue, res2 *pdb.Residue) float64 {
	minDist := distance(res1.Atoms[0], res2.Atoms[0])
	for _, a1 := range res1.Atoms {
		for _, a2 := range res2.Atoms {
			dist := distance(a1, a2)
			if dist < minDist {
				minDist = dist
			}
		}
	}

	return minDist
}

func maxResDist(pos int64, unp *uniprot.UniProt, p *pdb.PDB) (furthest residueDistance, err error) {
	fromRes := p.UniProtPositions[unp.ID][pos][0]

	for pos, residues := range p.UniProtPositions[unp.ID] {
		res := residues[0]
		dist := residuesDistance(fromRes, res)
		if dist > furthest.Distance {
			furthest = residueDistance{
				UnpPos:   pos,
				Residue:  res,
				Distance: dist,
			}
		}
	}

	return furthest, nil
}

func minGlycoDist(pos int64, unp *uniprot.UniProt, p *pdb.PDB) (closest residueDistance, closest2nd residueDistance, err error) {
	if len(unp.PTMs.Glycosilations) < 1 {
		return closest, closest2nd, errors.New("sequence does not have at least one glyco site")
	}

	// Glycosilation site residues (asparagines) in structure
	var glycoSites []residueDistance
	for _, site := range unp.PTMs.Glycosilations {
		if glycoSiteResidues, ok := p.UniProtPositions[unp.ID][site.Position]; ok {
			for _, res := range glycoSiteResidues {
				gs := residueDistance{
					UnpPos:   site.Position,
					Residue:  res,
					Distance: residuesDistance(res, p.UniProtPositions[unp.ID][pos][0]),
				}
				glycoSites = append(glycoSites, gs)
			}
		}
	}

	if len(glycoSites) < 1 {
		return closest, closest2nd, errors.New("crystal does not cover at least one glyco site")
	}

	// Sort sites by 3D distance
	sort.Slice(glycoSites, func(i, j int) bool {
		return glycoSites[i].Distance < glycoSites[j].Distance
	})

	closest = glycoSites[0]
	if len(glycoSites) > 1 {
		closest2nd = glycoSites[1]
	} else {
		closest2nd = residueDistance{
			UnpPos: 0, Residue: nil, Distance: math.NaN(),
		}
	}

	return closest, closest2nd, nil
}

func coveredGlycoSites(unp *uniprot.UniProt, p *pdb.PDB) (quantity int64) {
	for _, site := range unp.PTMs.Glycosilations {
		if _, ok := p.UniProtPositions[unp.ID][site.Position]; ok {
			quantity++
		}
	}

	return quantity
}
