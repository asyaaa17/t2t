// display_icarus.js


let selectedBlockElement = null;

// Инициализация отображения всех хромосом
function initAllChromosomesVisualization() {
    if (typeof references_by_id === 'undefined' || Object.keys(references_by_id).length === 0) {
        console.error("Нет данных для визуализации (references_by_id пуст).");
        document.getElementById("chart").innerHTML = "<div class='warning'>Нет данных для визуализации.</div>";
        return;
    }

    var chartContainer = document.getElementById("chart");

    for (var id in references_by_id) {
        var chrName = references_by_id[id];

        if (!chromosomes_len[chrName] || !contig_data[chrName]) {
            console.warn("Данные отсутствуют для хромосомы: " + chrName);
            continue;
        }

        var chrDiv = document.createElement("div");
        chrDiv.className = "chromosome_chart";
        chrDiv.id = "chart_" + chrName;
        chartContainer.appendChild(chrDiv);

        drawChromosome(chrName, chrDiv.id);
    }
}


// Функция для рисования каждой хромосомы с осью X
function drawChromosome(chrName, containerId) {
    var contigs = contig_data[chrName];
    var chrLength = chromosomes_len[chrName];

    var svgWidth = 1000;
    var svgHeight = 50;
    var marginLeft = 60;

    var svg = d3.select("#" + containerId)
        .append("svg")
        .attr("width", svgWidth)
        .attr("height", svgHeight + 30);  // Высота увеличена для оси X

    // Создание масштаба для оси X
    var xScale = d3.scale.linear()
        .domain([0, chrLength])
        .range([marginLeft, svgWidth - 40]);

    // Создание оси X
    var xAxis = d3.svg.axis()
        .scale(xScale)
        .orient('bottom')
        .ticks(6)
        .tickFormat(formatLongNumbers);  // 👈 используем форматирование Mbp/Kbp


    svg.append("defs").append("pattern")
        .attr("id", "centromerePattern")
        .attr("patternUnits", "userSpaceOnUse")
        .attr("width", 6)
        .attr("height", 6)
        .append("path")
        .attr("d", "M0,0 l6,6")
        .attr("stroke", "#888")  // серая штриховка
        .attr("stroke-width", 1);



    // Добавление оси X с тонким стилем
    svg.append('g')
        .attr('class', 'x axis')
        .attr('transform', 'translate(0,' + svgHeight + ')') // Перемещаем ось в нижнюю часть
        .call(xAxis)  // Применяем ось к графику
        .selectAll("path, line")
        .style("fill", "none")
        .style("stroke", "#000")  // Цвет линии оси
        .style("stroke-width", "0.5px");  // Устанавливаем тонкую ось


    var contigSetNames = Object.keys(contigs);
    var yOffset = 15;

    var rectGroup = svg.append("g").attr("class", "rectGroup");

    for (var i = 0; i < contigSetNames.length; i++) {
        var label = contigSetNames[i];
        var contigList = contigs[label];
        // enrich all contigs with objClass via createMiniItem
        for (var j = 0; j < contigList.length; j++) {
            createMiniItem(contigList[j], i, j, 0);  // i → curLane, j → numItem
        }

        rectGroup.selectAll("rect_" + label)
            .data(contigList)
            .enter()
            .append("rect")
            .each(function (d) {
                changeMisassembledStatus(d);
                if (d.name && d.name.toLowerCase().includes("centromere")) {
                    d.objClass = "centromere_block";
                }
            })
            .attr("x", function (d) { return xScale(d.start); })
            .attr("y", yOffset)
            .attr("width", function (d) { return Math.max(1, xScale(d.end) - xScale(d.start)); })
            .attr("height", 20)
            .attr("class", function (d) {
                if (d.misassembledEnds)
                    return "block end " + (d.objClass || "");
                if (!d.marks || d.contig_type)
                    return "block mainItem " + (d.objClass || "");
                return "block";
            })
            .attr("stroke", "#000")
            .attr("stroke-width", 0.5)
            .on("click", function (event, d) {  // D3 v5+ формат
                handleBlockClick(d, this);
            })
            .append("title")
            .text(function (d) {
                return d.name + ": " + d.start + " - " + d.end;
            });




        for (var j = 0; j < contigList.length; j++) {
            var block = contigList[j];  // <-- добавляем это

            createMiniItem(block, i, j, 0);  // i → curLane, j → numItem

            // рисуем прямоугольники
            rectGroup.append("rect")
                .datum(block)
                .each(function (d) {
                    changeMisassembledStatus(d);
                    if (d.name && d.name.toLowerCase().includes("centromere")) {
                        d.objClass = "centromere_block";
                    }
                })
                .attr("x", xScale(block.start))
                .attr("y", yOffset)
                .attr("width", Math.max(1, xScale(block.end) - xScale(block.start)))
                .attr("height", 20)
                .attr("class", function () {
                    if (block.misassembledEnds)
                        return "block end " + (block.objClass || "");
                    if (!block.marks || block.contig_type)
                        return "block mainItem " + (block.objClass || "");
                    return "block";
                })
                .attr("stroke", "#000")
                .attr("stroke-width", 0.5)
                .on("click", function () {
                    // вместо handleBlockClick(block, this);
                    // используй d3.select(this).datum() или d3.select(this).data()[0]
                    var blockData = d3.select(this).datum() || d3.select(this).data()[0];
                    handleBlockClick(blockData, this);
                })

                .append("title")
                .text(block.name + ": " + block.start + " - " + block.end);


        }




    }

    svg.append("text")
        .attr("x", 5)
        .attr("y", svgHeight / 2 + 5)
        .text(chrName)
        .attr("font-size", "14px")
        .attr("fill", "#333");

    svg.selectAll("rect")
        .sort(function (a, b) {
            return d3.ascending(a.start, b.start);
        });
}


function handleBlockClick(blockData, element) {
    const clicked = d3.select(element);

    // Если этот блок уже выделен — снять выделение
    if (selectedBlockElement === element) {
        clicked.attr("stroke-width", 0.5).attr("stroke", "#000");
        selectedBlockElement = null;
        clearInfoPanel();
        return;
    }

    // Снять выделение с предыдущего блока
    if (selectedBlockElement !== null) {
        d3.select(selectedBlockElement)
            .attr("stroke-width", 0.5)
            .attr("stroke", "#000");
    }

    // Выделить новый блок
    clicked.attr("stroke-width", 2.5).attr("stroke", "#000");
    selectedBlockElement = element;

    updateInfoPanel(blockData);
    console.log("blockData:", blockData);
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
        // Получаем длину хромосомы из глобального объекта chromosomes_len
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
            window.items.push(...blocks);   // <-- ВОТ ЭТА СТРОКА ГЛАВНАЯ!
        }
    }

    var itemsLayer = d3.select('body')
        .append('svg')
        .attr('id', 'items');
    window.itemsContainer = itemsLayer.append('g');

    // Теперь точно есть items → можно вызывать интерфейс
    if (typeof setupInterface !== 'undefined') {
        setupInterface();
    } else {
        console.warn("setupInterface() не найден");
    }

    console.log("Общее количество items:", items.length);
});