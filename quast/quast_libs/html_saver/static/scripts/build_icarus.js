console.log("[build_icarus.js] ✅ Script loaded");

var INTERLACE_BLOCKS_COLOR = false;
var INTERLACE_BLOCKS_VERT_OFFSET = false;
var isContigSizePlot = false;
var miniLanesHeight = 10;
var offsetsMiniY = [0, 1, 2];




function addAssemblyDescription(lanes) {
    for (var laneNum = 0; laneNum < lanes.length; laneNum++) {
        if (lanes[laneNum].label) {
            assemblyName = lanes[laneNum].label;
            var description = assemblyName + '\n';
            description += 'length: ' + assemblies_len[assemblyName] + '\n';
            description += 'contigs: ' + assemblies_contigs[assemblyName] + '\n';
            if (!isContigSizePlot)
                description += 'mis: ' + assemblies_misassemblies[assemblyName];
            else
                description += 'N50: ' + assemblies_n50[assemblyName];
            lanes[laneNum].description = description;
            if (!isContigSizePlot)
                lanes[laneNum].link = assemblies_links[assemblyName];
        }
    }
    return lanes;
}



function collapseLanes(chart) {
    console.log("[collapseLanes] called", chart);

    var lanes = [], items = [], laneId = 0, itemId = 0, groupId = 0;

    for (var assemblyName in chart.assemblies) {
        var lane = chart.assemblies[assemblyName];
        var currentLen = 0;
        var numItems = 0;
        var lastPosInLanes = [];
        var lastBlockSizesInLanes = [];
        var laneItems = [];

        function parseItem(block, fullInfo, misassembly) {
            block.misassembledEnds = '';
            block.lane = laneId;
            block.id = itemId;
            block.groupId = groupId;
            block.assembly = assemblyName;

            if (isContigSizePlot) {
                if (!fullInfo) {
                    block.corr_start = currentLen;
                    currentLen += block.size;
                    block.corr_end = currentLen;
                    block.fullContig = true;
                } else {
                    block.start_in_ref = block.corr_start;
                    block.end_in_ref = block.corr_end;
                    var start_in_contig = Math.min(block.start_in_contig, block.end_in_contig);
                    var end_in_contig = Math.max(block.start_in_contig, block.end_in_contig);
                    block.corr_start = currentLen + start_in_contig - 1;
                    block.corr_end = currentLen + end_in_contig;
                    block.notActive = true;
                    block.contig_type = fullInfo.contig_type;
                    block.mstype = misassembly ? misassembly.mstype : null;
                }
            }

            var nonOverlappingLaneId = 0;
            var blockSize = block.corr_end - block.corr_start;
            var minOverlap = Math.min(500, blockSize * 0.1);
            if (!isContigSizePlot) {
                for (nonOverlappingLaneId = 0; nonOverlappingLaneId < lastPosInLanes.length; nonOverlappingLaneId++) {
                    if (lastPosInLanes[nonOverlappingLaneId] - block.corr_start < Math.min(minOverlap, lastBlockSizesInLanes[nonOverlappingLaneId] * 0.1)) {
                        break;
                    }
                }
                block.nonOverlappingLane = nonOverlappingLaneId;
                if (nonOverlappingLaneId >= lastPosInLanes.length) {
                    lastPosInLanes.push(block.corr_end);
                    lastBlockSizesInLanes.push(blockSize);
                } else {
                    if (block.corr_end > lastPosInLanes[nonOverlappingLaneId]) {
                        lastPosInLanes[nonOverlappingLaneId] = block.corr_end;
                        lastBlockSizesInLanes[nonOverlappingLaneId] = blockSize;
                    }
                }
            }

            block.triangles = [];

            if (block.mis_ends) {
                var misassembled_ends = block.mis_ends.split(';');
                for (var num = 0; num < misassembled_ends.length; num++) {
                    var triangleItem = {};
                    triangleItem.name = block.name;
                    triangleItem.corr_start = block.corr_start;
                    triangleItem.corr_end = block.corr_end;
                    triangleItem.assembly = block.assembly;
                    triangleItem.id = itemId + num + 1;
                    triangleItem.lane = laneId;
                    triangleItem.nonOverlappingLane = block.nonOverlappingLane;
                    triangleItem.groupId = groupId;
                    triangleItem.misassembled = "True";
                    triangleItem.misassembledEnds = misassembled_ends[num];
                    triangleItem.misassemblies = block.misassemblies ? block.misassemblies.split(';')[num] : "";
                    triangleItem.objClass = 'misassembled';

                    block.triangles.push(triangleItem);
                    items.push(triangleItem);  //  ВАЖНО! Добавляем в items для отрисовки
                }
                itemId += misassembled_ends.length;
                numItems += misassembled_ends.length;
            }
            console.log(`[parseItem] block=${block.name || block.id} triangles=${block.triangles.length}`);

            itemId++;
            numItems++;
            return block;
        }

        for (var i = 0; i < lane.length; i++) {
            var block = lane[i];
            var newItem = parseItem(block);
            laneItems.push(newItem);
            groupId++;
        }

        for (var i = 0; i < laneItems.length; i++) {
            var item = laneItems[i];
            item.nonOverlappingLane = lastPosInLanes.length - item.nonOverlappingLane - 1;
            if (item.triangles) {
                for (var j = 0; j < item.triangles.length; j++) {
                    item.triangles[j].nonOverlappingLane = lastPosInLanes.length - item.triangles[j].nonOverlappingLane - 1;
                }
            }
            items.push(item);
        }

        lanes.push({
            id: laneId,
            label: assemblyName,
            maxLines: lastPosInLanes.length,
            isExpanded: false,
        });
        laneId++;
    }

    addAssemblyDescription(lanes);
    return { lanes: lanes, items: items };
}

function getMiniItems(items) {
    var result = [];
    var curLane = 0;
    var numItem = 0;
    var countSupplementary = 0;

    for (var i = 0; i < items.length; i++) {
        let block = items[i];
        if (block.lane != curLane) {
            numItem = 0;
            countSupplementary = 0;
        }

        result.push(createMiniItem(block, curLane, numItem, countSupplementary));
        curLane = block.lane;

        if (!block.notActive) numItem++;

        if (block.triangles && block.triangles.length > 0) {
            for (var j = 0; j < block.triangles.length; j++) {
                result.push(createMiniItem(block.triangles[j], curLane, numItem, countSupplementary));
                numItem++;
                countSupplementary++;
            }
        }
    }
    return result;
}

function createMiniItem(block, curLane, numItem, countSupplementary) {
    block.misassembled = block.misassemblies ? "True" : "False";
    let c = (block.misassembled == "True" ? "misassembled" : "");

    if (block.similar === "True") c += " similar";
    if (block.more_unaligned === "True" && block.misassembled === "True") {
        c = "mis_unaligned";
        block.misassemblies = "";
    } else if (block.more_unaligned === "True") {
        c = "correct_unaligned";
    }

    if (block.ambiguous === "True") c = "ambiguous";
    if (block.is_best !== "True" && block.ambiguous_alignments && block.ambiguous_alignments.length > 0) {
        c = "alternative";
    }

    if (INTERLACE_BLOCKS_COLOR) {
        c += ((numItem - countSupplementary) % 2 == 0 ? " odd" : "");
    }

    let text = '';
    if (block.marks) {
        text = block.marks;
        let marks = block.marks.split(', ');
        for (let m = 0; m < marks.length; m++) {
            c += " " + marks[m].toLowerCase();
        }
    }

    block.objClass = c;
    block.order = numItem - countSupplementary;

    if (block.triangles && block.triangles.length > 0) {
        block.misassembledEnds = block.triangles;
    }

    // Добавляем path (треугольник) для треугольных блоков
    if (block.misassembledEnds === "L") {
        let x = block.corr_start;
        let y = block.nonOverlappingLane * miniLanesHeight;
        let height = miniLanesHeight;

        block.path = `M${x},${y} L${x},${y + height} L${x - height / 1.5},${y + height / 2} Z`;
        console.log(`[createMiniItem] ↪️ Left triangle path=${block.path}`);
    } else if (block.misassembledEnds === "R") {
        let x = block.corr_end;
        let y = block.nonOverlappingLane * miniLanesHeight;
        let height = miniLanesHeight;

        block.path = `M${x},${y} L${x},${y + height} L${x + height / 1.5},${y + height / 2} Z`;
        console.log(`[createMiniItem] ↪️ Right triangle path=${block.path}`);
    } else if (block.misassembledEnds) {
        console.log(`[createMiniItem] ⚠️ Unknown misassembledEnds value: ${block.misassembledEnds}`);
    }



    //  console.log(`[createMiniItem] ${block.name || block.id || "?"} → objClass = ${block.objClass}, contig_type = ${block.contig_type}, misassemblies = ${block.misassemblies}, more_unaligned = ${block.more_unaligned}, similar = ${block.similar}, ambiguous = ${block.ambiguous}`);
}
