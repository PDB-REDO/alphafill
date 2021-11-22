import "core-js/stable";
import "regenerator-runtime/runtime";

import 'bootstrap';

import { Viewer } from "./molstar";

/* global page, PAGE_SIZE, STRUCTURE_COUNT */

class Pager {
	constructor(ul, pageCount, pageSize, tbl) {

		this.page = 1;
		this.lastPage = pageCount;
		this.pageSize = pageSize;
		this.tbl = tbl;

		this.ul = ul;
		this.buttons = [...ul.querySelectorAll('li')];
		this.pageButtonCount = this.buttons.length - 4;
		this.halfPageButtonCount = Math.round((this.pageButtonCount - 1) / 2);

		this.buttons.forEach(btn => {
			btn.addEventListener('click', () => {
				if (btn.getAttribute('data-page-id') !== undefined)
					this.turnPage(btn);
			});
		});

		this.updateButtons();
	}

	turnPage(btn) {
		let page = this.page;

		const pageID = btn.getAttribute('data-page-id');
		switch (pageID) {
			case 'previous':
				page = page - 1;
				break;
			case 'next':
				page = page + 1;
				break;
			default:
				page = parseInt(btn.textContent);
				break;
		}
		if (page > this.lastPage)
			page = this.lastPage;
		if (page < 1)
			page = 1;

		fetch(`./structure-table-page?page=${page}`, {
			credentials: "include"
		}).then(reply => {
			if (reply.ok)
				return reply.text();
			throw 'failed to fetch table';
		}).then(table => {
			this.tbl.innerHTML = table;
			this.page = page;
			this.updateButtons();
		}).catch(err => {
			console.log(err);
		});
	}

	updateButtons() {
		let pageForBtn1 = this.page - this.halfPageButtonCount;
		if (pageForBtn1 > this.lastPage - this.pageButtonCount)
			pageForBtn1 = this.lastPage - this.pageButtonCount;
		if (pageForBtn1 < 1)
			pageForBtn1 = 1;

		this.buttons[0].classList.toggle('disabled', this.page <= 1);
		this.buttons[1].style.display = pageForBtn1 > 1 ? '' : 'none';
	
		for (let i = 2; i < this.buttons.length - 2; ++i) {
			const pageNr = pageForBtn1 + i - 2;
			const btn = this.buttons[i];
			const a = btn.querySelector('a');
			a.textContent = pageNr;
			btn.classList.toggle('active', pageNr == this.page);
		}
	
		this.buttons[this.buttons.length - 2].style.display = pageForBtn1 < this.lastPage - this.pageButtonCount ? '' : 'none';
		this.buttons[this.buttons.length - 1].classList.toggle('disabled', this.page >= this.lastPage);
	}
}

window.addEventListener('load', () => {

	if (AF_ID != null) {
		const viewer = new Viewer(document.getElementById('app'), {
			layoutIsExpanded: false,
			layoutShowControls: false,
			viewportShowExpand: false,
			collapseLeftPanel: true,
			pixelScale: 1,
			pickScale: 0.25,
			pickPadding: 1,
		});
	
		viewer.loadStructureFromUrl(`/v1/aff/${AF_ID}`)
			.then(() => {
				const rows = document.querySelectorAll("tr.transplanted-row");
				[...rows].forEach(row => {
					row.addEventListener('click', () => {
						const asymID = row.getAttribute('data-asym-id');
						viewer.selectAsym(asymID);
					});
				});
			});
	}

	const rows = document.querySelectorAll('tr.structure-row');
	[...rows].forEach(row => {
		row.addEventListener('click', () => {
			const id = row.getAttribute('data-row-id');
			window.location = `./?id=${id}`;
		})
	});

	let lastPage = Math.round(STRUCTURE_COUNT / PAGE_SIZE);
	if (lastPage * PAGE_SIZE < STRUCTURE_COUNT)
		lastPage += 1;

	
	const ul = document.getElementById('stable-pager');
	const tbl = document.getElementById('structure-table');
	if (ul != null && tbl != null)
		new Pager(ul, lastPage, PAGE_SIZE, tbl);
})
