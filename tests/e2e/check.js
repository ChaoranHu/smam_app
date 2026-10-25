// End to end check of the app in a real browser, on a desktop and a phone screen.
// Usage: BASE_URL=http://127.0.0.1:8765 node tests/e2e/check.js [desktop|phone|both]
// Optional: SHOTS=dir for screenshots, CHECK_MATH=1 to check formula rendering (always on for https), MATHJAX_DIR=/usr/share/javascript/mathjax to serve
// MathJax locally when the CDN cannot be reached, FULL=1 to also run the slower checks.
const { chromium, devices } = require('playwright');
const fs = require('fs');
const path = require('path');

const BASE = process.env.BASE_URL || 'http://127.0.0.1:8765';
const SHOTS = process.env.SHOTS || 'e2e-shots';
const MJ = process.env.MATHJAX_DIR;
const FULL = process.env.FULL === '1';
const MATH = !!MJ || process.env.CHECK_MATH === '1' || /^https/.test(BASE);
const which = process.argv[2] || 'both';
fs.mkdirSync(SHOTS, { recursive: true });

const results = [];
function check(name, ok, detail = '') {
  results.push({ name, ok: !!ok, detail });
  console.log(`${ok ? 'PASS' : 'FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
}

async function idle(p, extra = 400) {
  await p.waitForTimeout(250);
  await p.waitForFunction(() => !document.documentElement.classList.contains('shiny-busy'), null, { timeout: 120000 });
  await p.waitForTimeout(extra);
}
const setSelect = (p, id, v) => p.evaluate(([id, v]) => {
  const el = document.getElementById(id);
  if (el.selectize) el.selectize.setValue(v); else { el.value = v; $(el).trigger('change'); }
}, [id, v]);
const setRadio = (p, id, v) => p.evaluate(([id, v]) => {
  const el = document.querySelector(`#${id} input[value="${v}"]`); el.click();
}, [id, v]);
const setNumber = (p, id, v) => p.evaluate(([id, v]) => { const el = document.getElementById(id); el.value = v; $(el).trigger('change'); }, [id, v]);
const text = (p, sel) => p.evaluate((sel) => { const e = document.querySelector(sel); return e ? e.innerText : ''; }, sel);
const errorsShown = (p) => p.evaluate(() => [...document.querySelectorAll('.shiny-output-error')]
  .filter(e => e.offsetParent !== null && !e.classList.contains('shiny-output-error-validation'))
  .map(e => e.id + ': ' + e.innerText.slice(0, 120)));
const noSideScroll = (p) => p.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth + 1);
const plotOk = (p, id) => p.evaluate((id) => { const i = document.querySelector(`#${id} img`); return !!i && i.naturalWidth > 50; }, id);

async function goTab(p, phone, label) {
  if (phone) {
    const open = await p.evaluate(() => !!document.querySelector('.navbar-collapse.show'));
    if (!open) { await p.click('nav button.navbar-toggler, nav button.navbar-toggle'); await p.waitForTimeout(500); }
  }
  await p.locator('a.nav-link', { hasText: label }).first().click();
  await idle(p, 700);
  if (phone) {
    const stillOpen = await p.evaluate(() => !!document.querySelector('.navbar-collapse.show'));
    check(`${label}: phone menu closes after picking the tab`, !stillOpen);
  }
}

async function mathRendered(p, sel) {
  return p.evaluate((sel) => {
    const root = document.querySelector(sel);
    if (!root) return { ok: false, why: 'missing' };
    const raw = /\$\$|\\\(|\\\)/.test(root.innerText);
    const n = root.querySelectorAll('.MathJax, .MathJax_Display, mjx-container, .MathJax_CHTML').length;
    return { ok: !raw && n > 0, why: `typeset elements ${n}, raw TeX left ${raw}` };
  }, sel);
}

async function run(kind) {
  const phone = kind === 'phone';
  const b = await chromium.launch(process.env.CHROMIUM ? { executablePath: process.env.CHROMIUM } : {});
  const ctx = await b.newContext(phone ? devices['iPhone 13'] : { viewport: { width: 1440, height: 950 }, acceptDownloads: true });
  if (MJ) {
    await ctx.route(/mathjax\.rstudio\.com\/latest\/(.*)/, (route) => {
      const rel = route.request().url().replace(/^.*\/latest\//, '').split('?')[0];
      route.fulfill({ path: path.join(MJ, rel) });
    });
  }
  const p = await ctx.newPage();
  const pageErrors = [];
  p.on('pageerror', e => pageErrors.push(e.message));
  const tag = (s) => `[${kind}] ${s}`;
  const shot = (n, full = true) => p.screenshot({ path: path.join(SHOTS, `${kind}_${n}.png`), fullPage: full });

  await p.goto(BASE, { waitUntil: 'load' });
  await p.waitForSelector('#data_boxes .bslib-value-box', { timeout: 90000 });
  await idle(p, 1500);

  // ---- Data ----
  check(tag('Data: default example loads with 4 value boxes'), (await p.locator('#data_boxes .bslib-value-box').count()) === 4);
  check(tag('Data: track plot draws'), await plotOk(p, 'plot_track'));
  check(tag('Data: no horizontal page scroll'), await noSideScroll(p));
  await shot('1_data');
  const examples = ['winter2011', 'summers', 'f109', 'simMR', 'summer2012'];
  for (const ex of examples) {
    await setSelect(p, 'example', ex); await idle(p, 1200);
    const boxes = await text(p, '#data_boxes');
    const errs = await errorsShown(p);
    check(tag(`Data: example ${ex} prepares`), /FIXES/i.test(boxes) && errs.length === 0, errs.join('; ') || boxes.replace(/\s+/g, ' ').slice(0, 80));
    if (ex === 'summers') check(tag('Data: summers example has segments'), /SEGMENTS\s*\n?\s*[1-9]/i.test(boxes), boxes.replace(/\s+/g, ' '));
    if (ex === 'f109') check(tag('Data: f109 season window builds segments'), /SEGMENTS\s*\n?\s*[1-9]/i.test(boxes), boxes.replace(/\s+/g, ' '));
  }
  // card tabs
  for (const [label, sel] of [['Coordinates over time', '#plot_coords img'], ['Sampling intervals', '#plot_steps img'],
                              ['Prepared data', '#table_prepared table'], ['Raw file', '#table_raw table'], ['Track', '#plot_track img']]) {
    await p.locator('.card-header-tabs a.nav-link', { hasText: label }).click(); await idle(p, 800);
    check(tag(`Data: "${label}" tab shows content`), (await p.locator(sel).count()) > 0);
  }
  // upload the CSV from the repository
  await setRadio(p, 'src', 'upload'); await idle(p);
  await p.setInputFiles('#file', path.join(__dirname, '..', '..', 'test.csv')); await idle(p, 1500);
  const up = await text(p, '#data_boxes');
  check(tag('Data: CSV upload maps columns and prepares'), /FIXES\s*\n?\s*(99|100)/i.test(up), up.replace(/\s+/g, ' ').slice(0, 60));
  check(tag('Data: upload picks time and coordinate columns'),
    (await p.inputValue('#time_col')) === 'time' &&
    (await p.evaluate(() => document.getElementById('coord_cols').selectize.getValue().join(','))) === 'east,north');
  // a bad mapping gives a message, not a crash
  await setSelect(p, 'coord_cols', ['time']); await idle(p);
  const msg = await text(p, '#data_boxes');
  check(tag('Data: invalid columns show a friendly message'), /cannot also|coordinate/i.test(msg), msg.slice(0, 80));
  await setSelect(p, 'coord_cols', ['east', 'north']); await idle(p);
  if (!phone) {
    const [dl] = await Promise.all([p.waitForEvent('download'), p.click('#dl_prepared')]);
    const f = await dl.path(); const lines = fs.readFileSync(f, 'utf8').trim().split('\n');
    check(tag('Data: prepared data downloads as CSV'), lines[0].includes('"time"') && lines.length > 50, `${lines.length} lines`);
  }
  await setRadio(p, 'src', 'example'); await idle(p);
  await setSelect(p, 'example', 'simMR'); await idle(p, 1200);

  // ---- Fit ----
  await goTab(p, phone, 'Fit model');
  check(tag('Fit: empty state before fitting'), /Choose a model and press Fit model/.test(await text(p, '#fit_results')));
  for (const [id, name] of [['MR', 'Moving-Resting (MR)'], ['MRME', 'Measurement Error (MRME)'], ['MRH', 'Handling (MRH)'],
                            ['MM', 'Moving-Moving (MM)'], ['BMME', 'Brownian Motion with Measurement Error']]) {
    await setSelect(p, 'model', id); await idle(p, 1200);
    const title = await text(p, '#model_math_title');
    check(tag(`Fit: model panel follows ${id}`), title.includes(name), title);
    if (MATH) {
      await p.waitForTimeout(1500);
      const m = await mathRendered(p, '#model_math');
      check(tag(`Fit: ${id} formulas render`), m.ok, m.why);
    }
    const nStart = await p.locator(`input[id^="start_${id}_"]:visible`).count();
    check(tag(`Fit: ${id} shows its starting values`), nStart === { MR: 3, MRME: 4, MRH: 5, MM: 4, BMME: 2 }[id], `${nStart} inputs`);
  }
  if (!phone) {
    await setSelect(p, 'model', 'MRME'); await idle(p, 1000);
    await shot('2_fit_model_mrme', false);
  }
  const fitModel = async (id, extra) => {
    await setSelect(p, 'model', id); await idle(p);
    if (extra) await extra();
    await p.click('#fit_go');
    await p.waitForFunction((id) => {
      const h = document.querySelector('#fit_results .card-header');
      const failed = /failed/.test((document.querySelector('#shiny-notification-panel') || {}).innerText || '');
      return failed || (h && h.innerText.includes('(' + id + ')') && !document.querySelector('#status_fit .status-card') &&
        document.querySelector('#est_table table'));
    }, id, { timeout: 240000 });
    await idle(p, 1000);
    const note = await text(p, '#shiny-notification-panel');
    if (/failed/.test(note)) check(tag(`Fit: ${id} runs`), false, note.slice(0, 200));
  };
  await setSelect(p, 'model', 'MRME'); await idle(p);
  await setRadio(p, 'MRME_variant', 'approx'); await idle(p);
  // versions deployed before the warning was added have no "very slow" label
  if (/very slow/.test(await text(p, '#MRME_variant'))) {
    check(tag('Fit: the slow approximation estimator shows a warning'), await p.locator('text=This estimator is very slow').isVisible());
  }
  await setRadio(p, 'MRME_variant', 'exact'); await idle(p);
  const fits = [['MR', null], ['BMME', null]];
  if (FULL) fits.push(['MM', null]);
  for (const [id, extra] of fits) {
    await fitModel(id, extra);
    const rows = await p.locator('#est_table tbody tr').count();
    const want = { MR: 3, MM: 4, BMME: 2 }[id];
    const hdr = await text(p, '#fit_results .card-header');
    check(tag(`Fit: ${id} gives ${want} estimates`), rows === want, `${rows} rows; ${hdr.replace(/\s+/g, ' ')}`);
    check(tag(`Fit: ${id} shows derived quantities and R code`), (await p.locator('#derived_table tbody tr').count()) > 0 &&
      /library\(smam\)/.test(await text(p, '#fit_results pre')));
    check(tag(`Fit: ${id} no error output`), (await errorsShown(p)).length === 0, (await errorsShown(p)).join('; '));
  }
  // MR with Hessian standard errors
  await fitModel('MR', () => setRadio(p, 'MR_likelihood', 'full'));
  const est = await text(p, '#est_table');
  check(tag('Fit: MR full likelihood shows SE and 95% CI'), /SE/.test(est) && /95% CI/.test(est));
  check(tag('Fit: convergence badge'), /converged/i.test(await text(p, '#fit_results .card-header')));
  await shot('3_fit_results');
  // switching the model warns that the estimates are from another model
  await setSelect(p, 'model', 'BMME'); await idle(p);
  check(tag('Fit: warning when the selected model differs from the fit'), /from the last fit/.test(await text(p, '#fit_mismatch')));
  await setSelect(p, 'model', 'MR'); await idle(p);
  check(tag('Fit: warning clears when the models match'), (await text(p, '#fit_mismatch')).trim() === '');
  if (!phone) {
    const [dl] = await Promise.all([p.waitForEvent('download'), p.click('#dl_estimates')]);
    const t = fs.readFileSync(await dl.path(), 'utf8');
    check(tag('Fit: estimates download'), /Estimate/.test(t) && /rate of leaving/.test(t) && !/&lambda;|<sub>/.test(t), t.split('\n')[1]);
  }
  // bootstrap SE for MRME approx is not offered; MM bootstrap is (FULL only)
  // cancel a long fit
  await setSelect(p, 'model', 'MRME'); await idle(p);
  await setRadio(p, 'MRME_variant', 'exact'); await idle(p);
  await p.click('#fit_go'); await p.waitForTimeout(2500);
  const running = await text(p, '#status_fit');
  check(tag('Fit: a running fit shows its status'), /Fitting/.test(running), running.slice(0, 60));
  if (!phone) await shot('4_fit_running', false);
  await p.click('#fit_cancel'); await idle(p, 1500);
  check(tag('Fit: cancel stops the fit'), (await text(p, '#status_fit')).trim() === '');
  check(tag('Fit: cancel does not report a failure'), !/failed/.test(await text(p, '#shiny-notification-panel')));

  // ---- States ----
  await goTab(p, phone, 'States');
  check(tag('States: uses the latest MR fit'), /lamM/.test(await text(p, '#state_fit_note')), await text(p, '#state_fit_note'));
  for (const method of ['fb', 'viterbi']) {
    await setRadio(p, 'state_method', method); await idle(p);
    await p.click('#state_go');
    await p.waitForFunction(() => document.querySelector('#state_track img') && !document.querySelector('#status_state .status-card'), null, { timeout: 180000 });
    await idle(p, 1500);
    const head = await text(p, '#state_results .card-header');
    check(tag(`States: MR ${method} draws track and path`), (await plotOk(p, 'state_track')) && (await plotOk(p, 'state_path')), head);
    if (method === 'fb') { check(tag('States: probability plot'), await plotOk(p, 'state_prob')); await shot('5_states'); }
    check(tag(`States: MR ${method} table`), (await p.locator('#state_table tbody tr').count()) > 0);
  }
  await setRadio(p, 'state_source', 'manual'); await idle(p);
  check(tag('States: manual parameters appear'), (await p.locator('input[id^="theta_MR_"]:visible').count()) === 3);
  await setRadio(p, 'state_source', 'fit'); await idle(p);

  // ---- Simulate ----
  await goTab(p, phone, 'Simulate');
  check(tag('Simulate: empty state'), /press Simulate/.test(await text(p, '#sim_results')));
  for (const id of ['MR', 'MRME', 'MRH', 'MM', 'BMME']) {
    await setSelect(p, 'sim_model', id); await idle(p);
    await p.click('#sim_go'); await idle(p, 1500);
    const ok = (await plotOk(p, 'sim_track')) && (await plotOk(p, 'sim_coords'));
    check(tag(`Simulate: ${id} track`), ok && (await errorsShown(p)).length === 0, await text(p, '#sim_results .card-header'));
  }
  await setNumber(p, 'sim_dim', 1); await p.click('#sim_go'); await idle(p, 1500);
  check(tag('Simulate: one dimension'), await plotOk(p, 'sim_track'));
  await setNumber(p, 'sim_dim', 3); await setRadio(p, 'sim_spacing', 'random'); await p.click('#sim_go'); await idle(p, 1500);
  check(tag('Simulate: three dimensions, irregular times'), (await plotOk(p, 'sim_track')) && /X3/.test(await text(p, '#sim_table')));
  await setNumber(p, 'sim_n', 0); await p.click('#sim_go'); await idle(p, 800);
  check(tag('Simulate: bad settings give a message'), /Check the simulation settings/.test(await text(p, '#shiny-notification-panel')));
  await setNumber(p, 'sim_n', 200); await setNumber(p, 'sim_dim', 2); await setSelect(p, 'sim_model', 'MR'); await idle(p);
  await p.click('#sim_go'); await idle(p, 1500);
  await shot('6_simulate');
  await p.click('#sim_use'); await idle(p, 2000);
  const active = await p.evaluate(() => document.querySelector('.navbar .nav-link.active').innerText);
  check(tag('Simulate: "Use this track" opens the Data tab with it'), /Data/.test(active) && (await p.inputValue('#example')) === 'simulated' && /FIXES\s*\n?\s*200/i.test(await text(p, '#data_boxes')));

  // ---- Guide ----
  await goTab(p, phone, 'Guide');
  check(tag('Guide: references listed'), (await p.locator('.guide ol li').count()) >= 9);
  if (MATH) {
    await p.waitForTimeout(2500);
    const m = await mathRendered(p, '.guide');
    check(tag('Guide: formulas render'), m.ok, m.why);
  }
  await shot('7_guide', false);
  for (const label of ['Data', 'Fit model', 'States', 'Simulate', 'Guide']) {
    await goTab(p, phone, label);
    check(tag(`${label}: no horizontal page scroll`), await noSideScroll(p));
  }
  check(tag('no JavaScript errors'), pageErrors.length === 0, pageErrors.join(' | ').slice(0, 200));
  await b.close();
}

(async () => {
  for (const k of (which === 'both' ? ['desktop', 'phone'] : [which])) {
    try { await run(k); } catch (e) { check(`[${k}] run finished`, false, e.message.split('\n')[0]); }
  }
  const failed = results.filter(r => !r.ok);
  console.log(`\n${results.length - failed.length} passed, ${failed.length} failed`);
  fs.writeFileSync(path.join(SHOTS, 'results.json'), JSON.stringify(results, null, 2));
  process.exit(failed.length ? 1 : 0);
})();
