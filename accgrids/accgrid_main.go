package main

import (
	"bufio"
	"compress/gzip"
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"gonum.org/v1/gonum/stat"
)

// aggrigate ascii grids over 30 years(aggRange) per 1 year (aggStep)
// Number of grids (N)
// generate N/aggStep - aggRange files

var crop = map[string]string{
	"WW_rcp85":      "wheatwinterwheat",
	"SM_rcp85":      "maizesilagemaize",
	"WRa_rcp85":     "rapewinterrape",
	"WW_rcp26":      "wheatwinterwheat",
	"SM_rcp26":      "maizesilagemaize",
	"WRa_rcp26":     "rapewinterrape",
	"WW_rcp85_420":  "wheatwinterwheat",
	"WRa_rcp85_420": "rapewinterrape",
	"SM_rcp85_420":  "maizesilagemaize",
	"WG_rcp85":      "barleywinterbarley",
	"WR_rcp85":      "ryewinterrye",
	"WW_rcp85_0":    "wheatwinterwheat",
	"WW_rcp85_1":    "wheatwinterwheat",
	"WW_rcp85_2":    "wheatwinterwheat",
	"WW_rcp85_3":    "wheatwinterwheat",
	"WW_rcp85_4":    "wheatwinterwheat",
	"WW_rcp85_5":    "wheatwinterwheat",
	"WG_rcp85_0":    "barleywinterbarley",
	"WG_rcp85_1":    "barleywinterbarley",
	"WG_rcp85_2":    "barleywinterbarley",
	"WG_rcp85_3":    "barleywinterbarley",
	"WG_rcp85_4":    "barleywinterbarley",
	"WG_rcp85_5":    "barleywinterbarley",
	"SM_rcp85_0":    "maizesilagemaize",
	"SM_rcp85_1":    "maizesilagemaize",
	"SM_rcp85_2":    "maizesilagemaize",
	"SM_rcp85_3":    "maizesilagemaize",
	"SM_rcp85_4":    "maizesilagemaize",
	"SM_rcp85_5":    "maizesilagemaize",
	"WR_rcp85_0":    "ryewinterrye",
	"WR_rcp85_1":    "ryewinterrye",
	"WR_rcp85_2":    "ryewinterrye",
	"WR_rcp85_3":    "ryewinterrye",
	"WR_rcp85_4":    "ryewinterrye",
	"WR_rcp85_5":    "ryewinterrye",
}

// var cropMask = map[string]string{
// 	"WW":  "de_grid_Rainfed_Wheat_GrowingArea_EU.csv",
// 	"WRa": "de_grid_Rainfed_Maize_GrowingArea_EU.csv",
// 	"SM":  "de_grid_Rainfed_Rapeseed_GrowingArea_EU.csv",
// }

var cropMask = map[string]string{
	"WW":  "CTM_17-19_mask_1000m_25832.asc",
	"WRa": "CTM_17-19_mask_1000m_25832.asc",
	"SM":  "CTM_17-19_mask_1000m_25832.asc",
	"WG":  "CTM_17-19_mask_1000m_25832.asc",
	"WR":  "CTM_17-19_mask_1000m_25832.asc",
}

var scenario = map[string]string{
	"WW_rcp85":      "rcp85",
	"SM_rcp85":      "rcp85",
	"WRa_rcp85":     "rcp85",
	"WW_rcp26":      "rcp26",
	"SM_rcp26":      "rcp26",
	"WRa_rcp26":     "rcp26",
	"WW_rcp85_420":  "rcp85",
	"WRa_rcp85_420": "rcp85",
	"SM_rcp85_420":  "rcp85",
	"WR_rcp85":      "rcp85",
	"WG_rcp85":      "rcp85",
	"WW_rcp85_0":    "rcp85",
	"WW_rcp85_1":    "rcp85",
	"WW_rcp85_2":    "rcp85",
	"WW_rcp85_3":    "rcp85",
	"WW_rcp85_4":    "rcp85",
	"WW_rcp85_5":    "rcp85",
	"WG_rcp85_0":    "rcp85",
	"WG_rcp85_1":    "rcp85",
	"WG_rcp85_2":    "rcp85",
	"WG_rcp85_3":    "rcp85",
	"WG_rcp85_4":    "rcp85",
	"WG_rcp85_5":    "rcp85",
	"SM_rcp85_0":    "rcp85",
	"SM_rcp85_1":    "rcp85",
	"SM_rcp85_2":    "rcp85",
	"SM_rcp85_3":    "rcp85",
	"SM_rcp85_4":    "rcp85",
	"SM_rcp85_5":    "rcp85",
	"WR_rcp85_0":    "rcp85",
	"WR_rcp85_1":    "rcp85",
	"WR_rcp85_2":    "rcp85",
	"WR_rcp85_3":    "rcp85",
	"WR_rcp85_4":    "rcp85",
	"WR_rcp85_5":    "rcp85",
}

var co2 = map[string]string{
	"WW_rcp85":      "-",
	"WRa_rcp85":     "-",
	"SM_rcp85":      "-",
	"WW_rcp85_420":  "420",
	"WRa_rcp85_420": "420",
	"SM_rcp85_420":  "420",
	"WW_rcp26":      "-",
	"SM_rcp26":      "-",
	"WRa_rcp26":     "-",
	"WR_rcp85":      "-",
	"WG_rcp85":      "-",
	"WW_rcp85_0":    "-",
	"WW_rcp85_1":    "-",
	"WW_rcp85_2":    "-",
	"WW_rcp85_3":    "-",
	"WW_rcp85_4":    "-",
	"WW_rcp85_5":    "-",
	"WG_rcp85_0":    "-",
	"WG_rcp85_1":    "-",
	"WG_rcp85_2":    "-",
	"WG_rcp85_3":    "-",
	"WG_rcp85_4":    "-",
	"WG_rcp85_5":    "-",
	"SM_rcp85_0":    "-",
	"SM_rcp85_1":    "-",
	"SM_rcp85_2":    "-",
	"SM_rcp85_3":    "-",
	"SM_rcp85_4":    "-",
	"SM_rcp85_5":    "-",
	"WR_rcp85_0":    "-",
	"WR_rcp85_1":    "-",
	"WR_rcp85_2":    "-",
	"WR_rcp85_3":    "-",
	"WR_rcp85_4":    "-",
	"WR_rcp85_5":    "-",
}

// hard coded setup map
// var setups = map[string][]int{
// 	"WW_rcp85":  {1, 41, 65, 81, 121, 129},
// 	"SM_rcp85":  {7, 47, 71, 87, 127, 135},
// 	"WRa_rcp85": {4, 44, 68, 84, 124, 132},
// 	"WW_rcp26":  {9, 17, 57, 89, 97},
// 	"SM_rcp26":  {15, 23, 63, 95, 103},
// 	"WRa_rcp26": {12, 20, 60, 92, 100},
// }

var setups = map[string][]int{
	"WW_rcp85_0": {137},
	"WW_rcp85_1": {157},
	"WW_rcp85_2": {169},
	"WW_rcp85_3": {177},
	"WW_rcp85_4": {197},
	"WW_rcp85_5": {201},
	"WG_rcp85_0": {138},
	"WG_rcp85_1": {158},
	"WG_rcp85_2": {170},
	"WG_rcp85_3": {178},
	"WG_rcp85_4": {198},
	"WG_rcp85_5": {202},
	"SM_rcp85_0": {140},
	"SM_rcp85_1": {160},
	"SM_rcp85_2": {172},
	"SM_rcp85_3": {180},
	"SM_rcp85_4": {200},
	"SM_rcp85_5": {204},
	"WR_rcp85_0": {139},
	"WR_rcp85_1": {159},
	"WR_rcp85_2": {171},
	"WR_rcp85_3": {179},
	"WR_rcp85_4": {199},
	"WR_rcp85_5": {203},
}

// da clim sims
// var setups = map[string][]int{
// 	"WW_rcp85":      {178, 182, 186, 190, 194, 198},
// 	"WRa_rcp85":     {180, 184, 188, 192, 196, 200},
// 	"SM_rcp85":      {181, 185, 189, 193, 197, 201},
// 	"WW_rcp85_420":  {202, 206, 210, 214, 218, 222},
// 	"WRa_rcp85_420": {204, 208, 212, 216, 220, 224},
// 	"SM_rcp85_420":  {205, 209, 213, 217, 221, 225},
// }

const inputFileformat = "%s_Yield_%d_%d.asc"

const inputFileformatPotET = "%s_Pot_ET_%d_%d.asc"
const inputFileformatPrecipSum = "%s_precip-sum_%d_%d.asc"
const inputFileformatTavgAvg = "%s_tavg-avg_%d_%d.asc"
const inputFileformatTmaxAvg = "%s_tmax-avg_%d_%d.asc"
const inputFileformatTminAvg = "%s_tmin-avg_%d_%d.asc"

const outStdfileTemplate = "avgYieldStd_%s_%s_%s_%d_%d.asc" // crop[setupId], scenario[setupId], co2, imageYear, index

// crop[setupId], scenario[setupId], co2, imageYear, index
var outDiffFileTemplate = map[refGridName]string{
	yields:    "avgYieldDiff_%s_%s_%s_%d_%d.asc",
	potET:     "avgPot_ETDiff_%s_%s_%s_%d_%d.asc",
	precipSum: "avgprecip-sumDiff_%s_%s_%s_%d_%d.asc",
	tavgAvg:   "avgtavg-avgDiff_%s_%s_%s_%d_%d.asc",
	tmaxAvg:   "avgtmax-avgDiff_%s_%s_%s_%d_%d.asc",
	tminAvg:   "avgtmin-avgDiff_%s_%s_%s_%d_%d.asc",
}

// crop[setupId], scenario[setupId], co2, imageYear, index
var outfileTemplate = map[refGridName]string{
	yields:    "avgYield_%s_%s_%s_%d_%d.asc",
	potET:     "Pot_ET_%s_%s_%s_%d_%d.asc",
	precipSum: "precip-sum_%s_%s_%s_%d_%d.asc",
	tavgAvg:   "tavg-avg_%s_%s_%s_%d_%d.asc",
	tmaxAvg:   "tmax-avg_%s_%s_%s_%d_%d.asc",
	tminAvg:   "tmin-avg_%s_%s_%s_%d_%d.asc",
}

var inputFolder = "./test"
var outFolder = "./agg_out"
var outFolders = map[string]string{}
var projectFolder = "."
var numConcurrent = 1
var aggRange uint = 30
var aggStep uint = 1
var cropId = "WW"
var startYear uint = 1971 // inclusive
var endYear uint = 2099   // inclusive
var withClimate = false
var rainfedLookup map[GridCoord]bool
var nodata = -9999.0

func main() {
	inputFolderPtr := flag.String("in", inputFolder, "path to input")
	outFolderPtr := flag.String("out", outFolder, "path to output")
	projectFolderPtr := flag.String("project", projectFolder, "path to output")
	concurrentPtr := flag.Int("concurrent", numConcurrent, "max concurrent execution")
	aggRangePtr := flag.Uint("aggRange", aggRange, "avarage of n years (default 30)")
	aggStepPtr := flag.Uint("aggStep", aggStep, "year jumps (default 1)")
	startyearPtr := flag.Uint("start", startYear, "start year (inclusive)")
	endyearPtr := flag.Uint("end", endYear, "end year (inclusive)")
	cropIdPtr := flag.String("crop", cropId, "crop id")
	withClimPtr := flag.Bool("climate", withClimate, "with climate")

	flag.Parse()
	inputFolder = *inputFolderPtr
	outFolder = *outFolderPtr
	projectFolder = *projectFolderPtr
	numConcurrent = *concurrentPtr
	withClimate = *withClimPtr

	aggRange = *aggRangePtr
	aggStep = *aggStepPtr
	cropId = *cropIdPtr
	startYear = *startyearPtr
	endYear = *endyearPtr

	if aggRange%2 > 0 {
		log.Fatal("aggRange should be an even number")
	}
	setupIds := []string{}
	for key := range setups {
		if strings.HasPrefix(key, cropId) {
			setupIds = append(setupIds, key)
			outFolders[key] = outFolder + strings.Trim(strings.Join(strings.Fields(fmt.Sprint(setups[key])), "_"), "[]")
		}
	}
	if len(setupIds) == 0 {
		log.Fatal("no setups found for crop")
	}
	rainfedSource := filepath.Join(projectFolder, cropMask[cropId])
	rainfedLookup = getMaskGridLookup(rainfedSource)

	startIdx := startYear + aggRange/2
	endIdx := endYear - aggRange/2 + 1

	gMinMax := newMinMax()
	outChan := make(chan MinMax)
	fileNameChan := make(chan fileMatch)
	final := make(chan MinMax)
	terminate := make(chan bool)
	go filenameCollector(fileNameChan, final, terminate)

	currRuns := 0

	for _, setupId := range setupIds {
		imageYearIndex := 0
		refGrids := map[refGridName][][]float64{}

		for imageYear := startIdx; imageYear < endIdx; imageYear = imageYear + aggStep {
			imageYearIndex++

			if imageYearIndex == 1 {
				// read reference grid (historical data) for diff maps
				lMinMax := calcAvgGrid(refGrids, withClimate, setupId, imageYear, imageYearIndex, fileNameChan, nil)
				gMinMax.setMinMax(&lMinMax)
			} else {
				go calcAvgGrid(refGrids, withClimate, setupId, imageYear, imageYearIndex, fileNameChan, outChan)

				currRuns++
				if currRuns >= numConcurrent {
					for currRuns >= numConcurrent {
						mMM := <-outChan
						currRuns--
						gMinMax.setMinMax(&mMM)

					}
				}
			}
		}
		for currRuns > 0 {
			mMM := <-outChan
			currRuns--
			gMinMax.setMinMax(&mMM)
		}
	}

	// send final min max values
	final <- gMinMax
	// wait for termination of meta file writing
	<-terminate
}

type fileMatch struct {
	ref      refGridName
	filename string
}

func filenameCollector(in chan fileMatch, final chan MinMax, out chan bool) {
	diffFilenames := []string{}
	yieldFilenames := map[refGridName][]string{}

	for {
		select {
		case globalMinMax := <-final:
			// run finished create meta files
			createMeta(yieldFilenames, diffFilenames, globalMinMax)
			// send finish signal to terminate
			out <- true
			return
		case filenameMatch := <-in:

			if filenameMatch.ref != numRefGridName {
				if _, ok := yieldFilenames[filenameMatch.ref]; !ok {
					yieldFilenames[filenameMatch.ref] = make([]string, 0, 100)
				}
				yieldFilenames[filenameMatch.ref] = append(yieldFilenames[filenameMatch.ref], filenameMatch.filename)
			} else if strings.Contains(filenameMatch.filename, "Diff") {
				diffFilenames = append(diffFilenames, filenameMatch.filename)
			}
		}
	}
}

type refGridName int

const (
	yields refGridName = iota
	stdYields
	potET
	precipSum
	tavgAvg
	tmaxAvg
	tminAvg
	numRefGridName
)

func calcAvgGrid(refGrids map[refGridName][][]float64, withClimate bool, setupId string, imageYear uint, imageYearIndex int, fileNameChan chan fileMatch, outC chan MinMax) (mMinMax MinMax) {
	type gridRefs struct {
		grid    [][]float64
		stdGrid [][][]float64
	}
	currentYearGrid := gridRefs{}
	currPotETGRid := gridRefs{}
	currPrecipSumGRid := gridRefs{}
	currTavgAvgGrid := gridRefs{}
	currTmaxAvgGrid := gridRefs{}
	currTminAvgGrid := gridRefs{}

	aggRangeHalf := aggRange / 2
	yearCounter := 0
	var header map[string]float64
	var stdeviationGrid [][]float64
	for _, setup := range setups[setupId] {

		// load grid - to grid buffer
		for imageIdx := imageYear - aggRangeHalf; imageIdx < imageYear+aggRangeHalf; imageIdx++ {
			// read grid file
			index := imageIdx - startYear + 1

			readFile := func(inFileformat string, grids *gridRefs) {
				filepath := filepath.Join(inputFolder, strconv.Itoa(setup), fmt.Sprintf(inFileformat, crop[setupId], imageIdx, index))
				file, err := os.Open(filepath)
				if err != nil {
					log.Fatal(err)
				}
				defer file.Close()
				scanner := bufio.NewScanner(file)

				if header == nil {
					header = readHeader(scanner, false)
				} else {
					// skip first lines
					readHeader(scanner, true)
				}
				if len(grids.grid) == 0 {
					cols := int(header["ncols"])
					rows := int(header["nrows"])
					grids.grid = make([][]float64, rows)
					grids.stdGrid = make([][][]float64, rows)
					for row := 0; row < rows; row++ {
						grids.grid[row] = make([]float64, cols)
						grids.stdGrid[row] = make([][]float64, cols)
					}
				}
				currRow := 0
				for scanner.Scan() {
					// sum up grid cells
					fields := strings.Fields(scanner.Text())
					for i, field := range fields {
						val, err := strconv.ParseFloat(field, 32)
						if err != nil {
							log.Fatal(err)
						}
						if val-nodata < 0.001 {
							grids.grid[currRow][i] = nodata
						} else {
							grids.grid[currRow][i] = grids.grid[currRow][i] + val
							if grids.stdGrid[currRow][i] == nil {
								grids.stdGrid[currRow][i] = make([]float64, 0, aggRange)
							}
							grids.stdGrid[currRow][i] = append(grids.stdGrid[currRow][i], val)
						}
					}
					currRow++
				}
			}
			readFile(inputFileformat, &currentYearGrid)

			if withClimate {
				readFile(inputFileformatPotET, &currPotETGRid)
				readFile(inputFileformatPrecipSum, &currPrecipSumGRid)
				readFile(inputFileformatTavgAvg, &currTavgAvgGrid)
				readFile(inputFileformatTmaxAvg, &currTmaxAvgGrid)
				readFile(inputFileformatTminAvg, &currTminAvgGrid)
			}
			yearCounter++
		}
		if stdeviationGrid == nil {
			stdeviationGrid = make([][]float64, int(header["nrows"]))
		}

		for rowIdx, row := range currentYearGrid.stdGrid {
			if stdeviationGrid[rowIdx] == nil {
				stdeviationGrid[rowIdx] = make([]float64, len(currentYearGrid.stdGrid[rowIdx]))
			}
			for colIdx := range row {
				if currentYearGrid.stdGrid[rowIdx][colIdx] != nil {
					val := stat.StdDev(currentYearGrid.stdGrid[rowIdx][colIdx], nil)
					if !math.IsNaN(val) {
						stdeviationGrid[rowIdx][colIdx] = stdeviationGrid[rowIdx][colIdx] + val
					} else {
						stdeviationGrid[rowIdx][colIdx] = nodata
					}
				} else {
					stdeviationGrid[rowIdx][colIdx] = nodata
				}
			}
		}
		cleanStdDev := func(gridRef *gridRefs) {
			for rowIdx, row := range gridRef.stdGrid {
				for colIdx := range row {
					gridRef.stdGrid[rowIdx][colIdx] = nil
				}
			}
		}
		cleanStdDev(&currentYearGrid)
		cleanStdDev(&currPotETGRid)
		cleanStdDev(&currPrecipSumGRid)
		cleanStdDev(&currTavgAvgGrid)
		cleanStdDev(&currTmaxAvgGrid)
		cleanStdDev(&currTminAvgGrid)
	}
	mMinMax = newMinMax()
	// calc average
	for rowIdx, row := range currentYearGrid.grid {
		for colIdx, col := range row {
			if currentYearGrid.grid[rowIdx][colIdx] != nodata {
				currentYearGrid.grid[rowIdx][colIdx] = col / float64(yearCounter)
				if stdeviationGrid[rowIdx][colIdx] != nodata {
					stdeviationGrid[rowIdx][colIdx] = stdeviationGrid[rowIdx][colIdx] / float64(len(setups[setupId]))

					mMinMax.setMaxV(int(stdeviationGrid[rowIdx][colIdx]), stdYields)
					mMinMax.setMinV(int(stdeviationGrid[rowIdx][colIdx]), stdYields)
				}
				mMinMax.setMaxV(int(currentYearGrid.grid[rowIdx][colIdx]), yields)
				mMinMax.setMinV(int(currentYearGrid.grid[rowIdx][colIdx]), yields)
				if withClimate {
					currPotETGRid.grid[rowIdx][colIdx] = currPotETGRid.grid[rowIdx][colIdx] / float64(yearCounter)
					currPrecipSumGRid.grid[rowIdx][colIdx] = currPrecipSumGRid.grid[rowIdx][colIdx] / float64(yearCounter)
					currTavgAvgGrid.grid[rowIdx][colIdx] = currTavgAvgGrid.grid[rowIdx][colIdx] / float64(yearCounter)
					currTmaxAvgGrid.grid[rowIdx][colIdx] = currTmaxAvgGrid.grid[rowIdx][colIdx] / float64(yearCounter)
					currTminAvgGrid.grid[rowIdx][colIdx] = currTminAvgGrid.grid[rowIdx][colIdx] / float64(yearCounter)

					mMinMax.setMaxV(int(currPotETGRid.grid[rowIdx][colIdx]), potET)
					mMinMax.setMaxV(int(currPrecipSumGRid.grid[rowIdx][colIdx]), precipSum)
					mMinMax.setMaxV(int(currTavgAvgGrid.grid[rowIdx][colIdx]), tavgAvg)
					mMinMax.setMaxV(int(currTmaxAvgGrid.grid[rowIdx][colIdx]), tmaxAvg)
					mMinMax.setMaxV(int(currTminAvgGrid.grid[rowIdx][colIdx]), tminAvg)
					mMinMax.setMinV(int(currPotETGRid.grid[rowIdx][colIdx]), potET)
					mMinMax.setMinV(int(currPrecipSumGRid.grid[rowIdx][colIdx]), precipSum)
					mMinMax.setMinV(int(currTavgAvgGrid.grid[rowIdx][colIdx]), tavgAvg)
					mMinMax.setMinV(int(currTmaxAvgGrid.grid[rowIdx][colIdx]), tmaxAvg)
					mMinMax.setMinV(int(currTminAvgGrid.grid[rowIdx][colIdx]), tminAvg)
				}
			}
		}
	}
	writeGrid := func(name refGridName, round int, useLookup bool, grid [][]float64) {
		// save new grid
		outFileName := filepath.Join(outFolders[setupId], fmt.Sprintf(outfileTemplate[name], crop[setupId], scenario[setupId], co2[setupId], imageYear, imageYearIndex))
		fout := writeAGridHeader(outFileName, header)
		writeFloatRows(fout, round, useLookup, grid)
		fout.Close()
		fileNameChan <- fileMatch{
			ref:      name,
			filename: outFileName,
		}
	}
	writeGrid(yields, 0, true, currentYearGrid.grid)

	// no diff grid for climate files
	if withClimate {
		writeGrid(potET, 2, false, currPotETGRid.grid)
		writeGrid(precipSum, 0, false, currPrecipSumGRid.grid)
		writeGrid(tavgAvg, 0, false, currTavgAvgGrid.grid)
		writeGrid(tmaxAvg, 0, false, currTmaxAvgGrid.grid)
		writeGrid(tminAvg, 0, false, currTminAvgGrid.grid)
	}

	writeDiffGrid := func(name refGridName, useLookup bool, grid [][]float64) {
		if refGrid, ok := refGrids[name]; ok {
			diffgrid := createDiff(refGrid, grid, nodata)

			outDiffFileName := filepath.Join(outFolders[setupId], fmt.Sprintf(outDiffFileTemplate[name], crop[setupId], scenario[setupId], co2[setupId], imageYear, imageYearIndex))
			fout := writeAGridHeader(outDiffFileName, header)
			writeIntRows(fout, useLookup, diffgrid)
			fout.Close()
			fileNameChan <- fileMatch{
				ref:      numRefGridName,
				filename: outDiffFileName,
			}
		} else {
			refGrids[name] = grid
		}
	}
	writeDiffGrid(yields, true, currentYearGrid.grid)

	// if withClimate {
	// 	writeDiffGrid(potET, currPotETGRid.grid)
	// 	writeDiffGrid(precipSum, currPrecipSumGRid.grid)
	// 	writeDiffGrid(tavgAvg, currTavgAvgGrid.grid)
	// 	writeDiffGrid(tmaxAvg, currTmaxAvgGrid.grid)
	// 	writeDiffGrid(tminAvg, currTminAvgGrid.grid)
	// }
	// save std grid
	outstdFileName := filepath.Join(outFolders[setupId], fmt.Sprintf(outStdfileTemplate, crop[setupId], scenario[setupId], co2[setupId], imageYear, imageYearIndex))
	foutStd := writeAGridHeader(outstdFileName, header)
	writeFloatRows(foutStd, 0, true, stdeviationGrid)
	foutStd.Close()
	fileNameChan <- fileMatch{
		ref:      stdYields,
		filename: outstdFileName,
	}

	if outC != nil {
		outC <- mMinMax
	}

	return mMinMax
}

func createDiff(refGrid, currentYearGrid [][]float64, nodata float64) (diffgrid [][]int) {

	diffgrid = make([][]int, len(currentYearGrid))

	calcDiffValue := func(hist, future int) int {
		diffVal := 0
		// catch diff by 0
		if hist == 0 {
			if future > 0 {
				diffVal = 100
			}
		} else {
			diffVal = (future - hist) * 100 / hist
		}
		// cap at +/-100%
		if diffVal > 100 {
			diffVal = 101
		}
		return diffVal
	}

	for rowIdx, row := range currentYearGrid {
		diffgrid[rowIdx] = make([]int, len(row))
		for colIdx := range row {
			if currentYearGrid[rowIdx][colIdx] == nodata {
				diffgrid[rowIdx][colIdx] = int(nodata)
			} else {
				diffgrid[rowIdx][colIdx] = calcDiffValue(int(refGrid[rowIdx][colIdx]), int(currentYearGrid[rowIdx][colIdx]))
			}
		}
	}
	return diffgrid
}

// readHeader and return a map of values, or skip the header
func readHeader(scanner *bufio.Scanner, skip bool) map[string]float64 {
	currRow := 0
	if skip {
		for currRow < 6 && scanner.Scan() {
			currRow++
		}
		return nil
	}
	outMap := make(map[string]float64, 6)

	for currRow < 6 && scanner.Scan() {
		currRow++
		texts := strings.Fields(scanner.Text())
		val, err := strconv.ParseFloat(texts[1], 32)
		if err != nil {
			log.Fatal(err)
		}
		outMap[texts[0]] = val
	}
	return outMap
}

// Fout combined file writer
type Fout struct {
	file    *os.File
	gfile   *gzip.Writer
	fwriter *bufio.Writer
}

// Write string to zip file
func (f Fout) Write(s string) {
	f.fwriter.WriteString(s)
}

// Close file writer
func (f Fout) Close() {
	f.fwriter.Flush()
	// Close the gzip first.
	f.gfile.Close()
	f.file.Close()
}

func writeAGridHeader(name string, header map[string]float64) (fout Fout) {
	cornerX := header["xllcorner"]
	cornery := header["yllcorner"]
	novalue := int(nodata)
	cellsize := header["cellsize"]
	nCol := int(header["ncols"])
	nRow := int(header["nrows"])
	// create an ascii file, which contains the header
	makeDir(name)
	file, err := os.OpenFile(name+".gz", os.O_CREATE|os.O_TRUNC|os.O_WRONLY, 0600)
	if err != nil {
		log.Fatal(err)
	}

	gfile := gzip.NewWriter(file)
	fwriter := bufio.NewWriter(gfile)
	fout = Fout{file, gfile, fwriter}

	fout.Write(fmt.Sprintf("ncols %d\n", nCol))
	fout.Write(fmt.Sprintf("nrows %d\n", nRow))
	fout.Write(fmt.Sprintf("xllcorner     %f\n", cornerX))
	fout.Write(fmt.Sprintf("yllcorner     %f\n", cornery))
	fout.Write(fmt.Sprintf("cellsize      %f\n", cellsize))
	fout.Write(fmt.Sprintf("NODATA_value  %d\n", novalue))

	return fout
}

func writeFloatRows(fout Fout, round int, useIrrgLookup bool, grid [][]float64) {

	for irow, row := range grid {
		for icol, col := range row {
			val := col
			if val-(nodata) < 0.001 {
				fout.Write(strconv.Itoa(int(val)))
			} else {
				if useIrrgLookup {
					if _, ok := rainfedLookup[GridCoord{irow, icol}]; !ok {
						val = -1
					}
				}
				if round == 0 {
					fout.Write(strconv.Itoa(int(math.Round(val))))
				} else {
					fout.Write(strconv.FormatFloat(val, 'f', round, 64))
				}
			}
			fout.Write(" ")

		}
		fout.Write("\n")
	}
}

func writeIntRows(fout Fout, useIrrgLookup bool, grid [][]int) {

	for irow, row := range grid {
		for icol, col := range row {
			val := col
			if val != int(nodata) && useIrrgLookup {
				if _, ok := rainfedLookup[GridCoord{irow, icol}]; !ok {
					val = -101
				}
			}
			fout.Write(strconv.Itoa(val))

			fout.Write(" ")
		}
		fout.Write("\n")
	}
}

func makeDir(outPath string) {
	dir := filepath.Dir(outPath)
	if _, err := os.Stat(dir); os.IsNotExist(err) {
		if err := os.MkdirAll(dir, os.ModePerm); err != nil {
			log.Fatalf("ERROR: Failed to generate output path %s :%v", dir, err)
		}
	}
}

func createMeta(yieldFilelist map[refGridName][]string, diffFileList []string, globalMinMax MinMax) {

	metaMap := make(map[refGridName]metaSetup, numRefGridName)
	for i := yields; i < numRefGridName; i++ {
		if i == stdYields {
			metaMap[i] = newStdMetaSet(globalMinMax.minVal[i], globalMinMax.maxVal[i])
		} else if i == precipSum {
			metaMap[i] = newYieldMetaSetup(globalMinMax.minVal[i], globalMinMax.maxVal[i], "jet_r")
		} else if i == tavgAvg || i == tmaxAvg || i == tminAvg {
			metaMap[i] = newTempMetaSetup()
		} else if i == potET {
			mm := newYieldMetaSetup(globalMinMax.minVal[i], globalMinMax.maxVal[i], "coolwarm")
			mm.colorList = append(metaMap[i].colorList, "royalblue", "mediumslateblue", "crimson")
			mm.colorListType = "LinearSegmented"
			metaMap[i] = mm
		} else if i == yields {
			metaMap[i] = newYieldMetaSetup(globalMinMax.minVal[i], globalMinMax.maxVal[i], "YlOrBr")
		} else {
			metaMap[i] = newYieldMetaSetup(globalMinMax.minVal[i], globalMinMax.maxVal[i], "jet")
		}
	}

	diffMeta := newDiffMetaSet()
	for ref, filenames := range yieldFilelist {
		for _, filename := range filenames {
			writeMetaFile(filename, metaMap[ref])
		}
	}
	for _, filename := range diffFileList {
		writeMetaFile(filename, diffMeta)
	}
}

type metaSetup struct {
	colormap string
	maxValue int
	minValue int

	colorList     []string
	colorListType string
	minColor      string
	border        bool
}

func newDiffMetaSet() metaSetup {
	return metaSetup{
		colormap:      "RdYlGn",
		maxValue:      101,
		minValue:      -101,
		colorList:     []string{},
		colorListType: "",
		minColor:      "lightgrey",
		border:        false,
	}
}
func newTempMetaSetup() metaSetup {
	return metaSetup{
		colormap:      "temperature",
		maxValue:      56,
		minValue:      -46,
		colorList:     []string{},
		colorListType: "",
		minColor:      "lightgrey",
		border:        false,
	}
}
func newYieldMetaSetup(minValue, maxValue int, colorMap string) metaSetup {
	return metaSetup{
		colormap:      colorMap,
		maxValue:      maxValue,
		minValue:      minValue,
		colorList:     []string{},
		colorListType: "",
		minColor:      "lightgrey",
		border:        false,
	}
}
func newStdMetaSet(minValue, maxValue int) metaSetup {
	return metaSetup{
		colormap:      "cool",
		maxValue:      maxValue,
		minValue:      minValue,
		colorList:     []string{},
		colorListType: "",
		minColor:      "lightgrey",
		border:        false,
	}
}

func writeMetaFile(gridFilePath string, setup metaSetup) {
	metaFilePath := gridFilePath + ".meta"
	makeDir(metaFilePath)
	file, err := os.OpenFile(metaFilePath, os.O_CREATE|os.O_TRUNC|os.O_WRONLY, 0600)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	if setup.colormap != "" {
		file.WriteString(fmt.Sprintf("colormap: '%s'\n", setup.colormap))
	}
	file.WriteString(fmt.Sprintf("maxValue: %d\n", setup.maxValue))
	file.WriteString(fmt.Sprintf("minValue: %d\n", setup.minValue))
	if len(setup.minColor) > 0 {
		file.WriteString(fmt.Sprintf("minColor: %s\n", setup.minColor))
	}
	if len(setup.colorList) > 0 {
		file.WriteString("colorlist: \n")
		for _, item := range setup.colorList {
			file.WriteString(fmt.Sprintf(" - '%s'\n", item))
		}
	}
	if len(setup.colorListType) > 0 {
		file.WriteString(fmt.Sprintf("colorlisttype: %s\n", setup.colorListType))
	}
	if setup.border {
		file.WriteString("border: True\n")
		file.WriteString("removeEmptyColumns: False\n")
	}
}

type MinMax struct {
	minVal map[refGridName]int
	maxVal map[refGridName]int
}

func newMinMax() MinMax {

	mm := MinMax{
		minVal: make(map[refGridName]int, numRefGridName),
		maxVal: make(map[refGridName]int, numRefGridName),
	}
	for i := yields; i < numRefGridName; i++ {
		mm.minVal[i] = -1
		mm.maxVal[i] = -1
	}
	return mm
}

func (m *MinMax) setMinV(min int, fieldName refGridName) {
	curValue := m.minVal[fieldName]
	if curValue > min || curValue < 0 {
		m.minVal[fieldName] = min
	}
}

func (m *MinMax) setMaxV(max int, fieldName refGridName) {
	curValue := m.maxVal[fieldName]
	if curValue < max {
		m.maxVal[fieldName] = max
	}
}

func (m *MinMax) setMinMax(other *MinMax) {
	for i := yields; i < numRefGridName; i++ {
		m.setMinV(other.minVal[i], i)
		m.setMaxV(other.maxVal[i], i)
	}
}

// GridCoord tuple of positions
type GridCoord struct {
	row int
	col int
}

func getMaskGridLookup(gridsource string) map[GridCoord]bool {
	lookup := make(map[GridCoord]bool)

	sourcefile, err := os.Open(gridsource)
	if err != nil {
		log.Fatal(err)
	}
	defer sourcefile.Close()

	scanner := bufio.NewScanner(sourcefile)
	header := readHeader(scanner, false)
	nodataVal := int64(header["NODATA_value"])

	row := 0
	for scanner.Scan() {
		line := scanner.Text()
		tokens := strings.Fields(line)
		for col := 0; col < len(tokens); col++ {
			landUsed, _ := strconv.ParseInt(tokens[col], 10, 64)
			if landUsed > 0 && landUsed != nodataVal {
				lookup[GridCoord{int(row), int(col)}] = true
			}
		}
		row++
	}
	return lookup
}

// func getMaskGridLookupFromCSV(gridsource string) map[GridCoord]bool {
// 	lookup := make(map[GridCoord]bool)

// 	sourcefile, err := os.Open(gridsource)
// 	if err != nil {
// 		log.Fatal(err)
// 	}
// 	defer sourcefile.Close()
// 	firstLine := true
// 	colID := -1
// 	rowID := -1
// 	irrID := -1
// 	scanner := bufio.NewScanner(sourcefile)
// 	for scanner.Scan() {
// 		line := scanner.Text()
// 		tokens := strings.Split(line, ",")
// 		if firstLine {
// 			firstLine = false
// 			// Column,Row,latitude,longitude,irrigation
// 			for index, token := range tokens {
// 				if token == "Column" {
// 					colID = index
// 				}
// 				if token == "Row" {
// 					rowID = index
// 				}
// 				if token == "irrigation" {
// 					irrID = index
// 				}
// 			}
// 		} else {
// 			col, _ := strconv.ParseInt(tokens[colID], 10, 64)
// 			row, _ := strconv.ParseInt(tokens[rowID], 10, 64)
// 			irr, _ := strconv.ParseInt(tokens[irrID], 10, 64)
// 			if irr > 0 {
// 				lookup[GridCoord{int(row), int(col)}] = true
// 			}
// 		}
// 	}
// 	return lookup
// }
