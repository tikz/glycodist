package main

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"log"
	"os"
	"sort"
)

func init() {
	makeDirs()
}

func main() {
	fOut, err := os.Create("result.csv")
	if err != nil {
		log.Fatal(err)
	}
	defer fOut.Close()

	writer := csv.NewWriter(fOut)
	defer writer.Flush()
	writer.Write([]string{"UniProtID", "PDB ID", "Pos", "ClosestGlycoPos", "ClosestGlycoDistance", "2ndClosestGlycoPos", "2ndClosestGlycoDistance", "FurthestResPos", "FurthestResDistance"})

	f, _ := os.Open("unps.txt")
	scanner := bufio.NewScanner(f)

	for scanner.Scan() {
		unpID := scanner.Text()
		unp, err := loadUniProt(unpID)
		if err != nil {
			continue
		}

		for i, pdbID := range unp.PDBIDs {
			if unp.PDBIDsCoverage[i] > 0.9 && len(unp.PTMs.Glycosilations) > 0 {
				pdb, err := loadPDB(pdbID)
				if err != nil {
					continue
				}

				// Distance to closest N-glyco site for all positions
				// in chains of the PDB that correspond to the UniProt ID

				// Sort map keys
				var coveredPos []int
				for pos := range pdb.UniProtPositions[unp.ID] {
					coveredPos = append(coveredPos, int(pos))
				}
				sort.Ints(coveredPos)

				for _, pos := range coveredPos {
					site, site2nd, err := minGlycoDist(int64(pos), unp, pdb)
					if err != nil {
						log.Println(unp.ID, pdb.ID, err)
						continue
					}

					furthestRes, err := maxResDist(int64(pos), unp, pdb)
					if err != nil {
						log.Println(unp.ID, pdb.ID, err)
						continue
					}

					line := []string{
						unp.ID,
						pdb.ID,
						fmt.Sprintf("%d", pos),
						fmt.Sprintf("%d", site.UnpPos),
						fmt.Sprintf("%f", site.Distance),
						fmt.Sprintf("%d", site2nd.UnpPos),
						fmt.Sprintf("%f", site2nd.Distance),
						fmt.Sprintf("%d", furthestRes.UnpPos),
						fmt.Sprintf("%f", furthestRes.Distance),
					}
					writer.Write(line)
					fmt.Println(line)
				}
			}
		}
	}
}
