// display_icarus.js


let selectedBlockElement = null;
let floatingInfoBox = null;


function initAllChromosomesVisualization() {
    if (typeof references_by_id === 'undefined' || Object.keys(references_by_id).length === 0) {
        console.error("Нет данных для визуализации (references_by_id пуст).");
        document.getElementById("chart").innerHTML = "<div class='warning'>Нет данных для визуализации.</div>";
        return;
    }

    var chartContainer = document.getElementById("chart");

    for (var id in references_by_id) {
        var chrName = references_by_id[id];


        var chrDiv = document.createElement("div");
        chrDiv.className = "chromosome_chart";
        chrDiv.id = "chart_" + chrName;
        chartContainer.appendChild(chrDiv);

        drawChromosome(chrName, chrDiv.id);
    }
}


function drawChromosome(chrName, containerId) {
    var contigs = contig_data[chrName];
    var chrLength = chromosomes_len[chrName];

    var svgWidth = 1000;
    var laneHeight = 22; 
    var marginTop = 30; 
    var marginLeft = 60;

    var assemblyLabels = Object.keys(contigs);
    var svgHeight = laneHeight * assemblyLabels.length + marginTop + 20;

    var colorScale = d3.scale.category10();

    var svg = d3.select("#" + containerId)
        .append("svg")
        .attr("width", svgWidth)
        .attr("height", svgHeight);

    // pattern for centromere
    svg.append("defs").append("pattern")
        .attr("id", "centromerePattern")
        .attr("patternUnits", "userSpaceOnUse")
        .attr("width", 6)
        .attr("height", 6)
        .append("path")
        .attr("d", "M0,0 l6,6")
        .attr("stroke", "#888")
        .attr("stroke-width", 1);

    // X-axis scale
    var xScale = d3.scale.linear()
        .domain([0, chrLength])
        .range([marginLeft, svgWidth - 40]);

    // X-axis
    var xAxis = d3.svg.axis()
        .scale(xScale)
        .orient('bottom')
        .ticks(6)
        .tickFormat(formatLongNumbers);

    svg.append('g')
        .attr('class', 'x axis')
        .attr('transform', 'translate(0,' + (svgHeight - 20) + ')')
        .call(xAxis)
        .selectAll("path, line")
        .style("fill", "none")
        .style("stroke", "#000")
        .style("stroke-width", "0.5px");

    assemblyLabels.forEach(function (label, idx) {
        var contigList = contigs[label];

        // name assemblies 
        svg.append("text")
            .attr("x", 0)
            .attr("y", marginTop + laneHeight * idx + 15)
            .text(label)
            .attr("font-size", "12px")
            .attr("fill", colorScale(idx));

            svg.selectAll("rect_" + label)
            .data(contigList)
            .enter()
            .append("rect")


            .attr("x", function (d) { return xScale(d.start); })
            .attr("y", marginTop + laneHeight * idx)
            .attr("width", function (d) { return Math.max(1, xScale(d.end) - xScale(d.start)); })
            .attr("height", laneHeight - 4)
            .attr("fill", function (d) {
                // for centromere
                if (d.name && d.name.toLowerCase().includes("centromere")) {
                    return "url(#centromerePattern)";
                }
                return colorScale(idx);
            })

            .attr("class", function (d) {
                let laneClass = "lane_" + label;  
                let centromereClass = "";

                if (d.name && d.name.toLowerCase().includes("centromere")) {
                    centromereClass = " centromere_block";
                }

                if (d.misassembledEnds)
                    return "block end " + (d.objClass || "") + " " + laneClass + centromereClass;
                if (!d.marks || d.contig_type)
                    return "block mainItem " + (d.objClass || "") + " " + laneClass + centromereClass;
                return "block " + laneClass + centromereClass;
            })


            .attr("stroke", "#000")
            .attr("stroke-width", 0.5)
            .on("click", function (d, i) { 
                console.log("[Клик по rect]:", d);
                handleBlockClick(d, this);
            })


            .append("title")
            .text(function (d) {
                return d.name + ": " + d.start + " - " + d.end;
            });


        svg.selectAll("triangle_" + label)
            .data(contigList.filter(d => (d.misassembled === "True") && (d.mis_ends === "L" || d.mis_ends === "R")))

            .enter()
            .append("path")
            .attr("d", function (d) {
                console.log("xStart:", xScale(d.start), "xEnd:", xScale(d.end), "mis_ends:", d.mis_ends);

                const xStart = xScale(d.start);
                const xEnd = xScale(d.end);
                const yTop = marginTop + laneHeight * idx;
                const yMid = yTop + (laneHeight - 4) / 2;
                const size = 6;

                if (d.mis_ends === "L") {
                    return `M${xStart},${yMid - size} L${xStart},${yMid + size} L${xStart + size},${yMid} Z`;
                } else if (d.mis_ends === "R") {
                    return `M${xEnd},${yMid - size} L${xEnd},${yMid + size} L${xEnd - size},${yMid} Z`;
                }

                return null;
            })
            .attr("class", d => "end_triangle " + (d.objClass || ""))

            .attr("class", "mis_triangle")

    });


    svg.append("text")
        .attr("x", 5)
        .attr("y", 18)
        .text(chrName)
        .attr("font-size", "16px")
        .attr("fill", "#333");
}




function handleBlockClick(blockData, element) {
    const clicked = d3.select(element);
    
    if (selectedBlockElement === element) {
        clicked.attr("stroke-width", 0.5).attr("stroke", "#000");
        selectedBlockElement = null;
        clearInfoPanel();
        if (floatingInfoBox) {
            floatingInfoBox.remove();
            floatingInfoBox = null;
        }
        return;
    }

    if (selectedBlockElement !== null) {
        d3.select(selectedBlockElement)
            .attr("stroke-width", 0.5)
            .attr("stroke", "#000");
    }

    clicked.attr("stroke-width", 2.5).attr("stroke", "#000");
    selectedBlockElement = element;

    updateInfoPanel(blockData);
    console.log("blockData:", blockData);

    // dlete 
    if (floatingInfoBox) {
        floatingInfoBox.remove();
        floatingInfoBox = null;
    }

    const svg = d3.select(element.ownerSVGElement);
    const x = +clicked.attr("x");
    const y = +clicked.attr("y");
    const lines = [
        `${blockData.name || "Contig"}, Start: ${formatLongNumbers(blockData.start)}, End: ${formatLongNumbers(blockData.end)}`
    ];

    // dinamic
    const tempSvg = d3.select("body").append("svg").attr("visibility", "hidden");
    let maxWidth = 0;
    lines.forEach(line => {
        const tempText = tempSvg.append("text").attr("font-size", "12px").text(line);
        const length = tempText.node().getComputedTextLength();
        if (length > maxWidth) maxWidth = length;
        tempText.remove();
    });
    tempSvg.remove();

    const boxWidth = maxWidth + 20;
    const lineHeight = 14;
    const boxHeight = lines.length * lineHeight + 8;

    // check
    const svgNode = svg.node();
    const svgWidth = +svgNode.getAttribute("width");
    let adjustedX = x;
    if (x + boxWidth > svgWidth - 10) {
        adjustedX = svgWidth - boxWidth - 10;
    }

    // create box
    floatingInfoBox = svg.append("g")
        .attr("class", "floating-info");

    floatingInfoBox.append("rect")
        .attr("x", adjustedX)
        .attr("y", y - boxHeight - 5)
        .attr("width", boxWidth)
        .attr("height", boxHeight)
        .attr("fill", "#fff")
        .attr("stroke", "#333")
        .attr("stroke-width", 1)
        .attr("rx", 4)
        .attr("ry", 4);

    lines.forEach((line, i) => {
        floatingInfoBox.append("text")
            .attr("x", adjustedX + 6)
            .attr("y", y - boxHeight + 14 * i + 8)
            .attr("font-size", "12px")
            .attr("fill", "#000")
            .text(line);
    });
}


function formatLongNumbers(d) {
    if (d >= 1000000) {
        return (d / 1000000).toFixed(2).replace(/\.0+$/, '') + " Mbp";
    } else if (d >= 1000) {
        return (d / 1000).toFixed(1).replace(/\.0$/, '') + " Kbp";
    } else {
        return d + " bp";
    }
}

function updateInfoPanel(blockData) {
    var infoPanel = document.getElementById("block_info");
    if (infoPanel) {
        let chrSize = chromosomes_len[blockData.name];
        let sizeStr = chrSize ? " (" + formatLongNumbers(chrSize) + ")" : "";

        infoPanel.innerHTML =
            `<b>${blockData.name || "?"}${sizeStr}</b><br/>
            Start: ${blockData.start || "?"}<br/>
            End: ${blockData.end || "?"}<br/>
            Contig Type: ${blockData.contig_type || "Unknown"}<br/>
            Is Best: ${blockData.is_best || "False"}<br/>
            Misassemblies: ${blockData.misassemblies || "None"}`;
    } else {
        console.warn("info_panel not found in DOM");
    }
}

function clearInfoPanel() {
    var infoPanel = document.getElementById("block_info");
    if (infoPanel) {
        infoPanel.innerHTML = "";
    }
}


document.addEventListener("DOMContentLoaded", function () {
    initAllChromosomesVisualization();
    window.items = [];
    window.items = [];
    for (var chr in contig_data) {
        var contigSets = contig_data[chr];
        for (var label in contigSets) {
            var blocks = contigSets[label];
            window.items.push(...blocks);   

        }
    }

    var itemsLayer = d3.select('body')
        .append('svg')
        .attr('id', 'items');
    window.itemsContainer = itemsLayer.append('g');

    if (typeof setupInterface !== 'undefined') {
        setupInterface();
    } else {
        console.warn("setupInterface() не найден");
    }

    console.log("Общее количество items:", items.length);
});
