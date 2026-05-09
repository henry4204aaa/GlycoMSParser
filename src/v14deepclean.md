# V14_man_deepreview.md — 

**Audit date:** 2026-04-12 - 2026-04-21
**Auditor:** Henry
**Reference:**: `V14_FUNCTION_AUDIT.md` and ``
**Working on:** `mspfileloaderv14.py` (10,428 lines, current on-disk state post-Phase 2 cleanup)

@Henry is comment from Henry

### Blocks chained and organized, functions reviewed

# diagnostics block
## line 33 - 100
@Henry: This blocks allow --diag to test essential module installation status on Windows.
Still from user's report, recent version sometimes fails with a flashout so they can't catch errors.

# json type/schema enforcement (GlycoMSP v1)
## line 106 - 220
@Henry: This block is a fixed, force applying rules to align 3 type of jsons with proper headers and fields.
Mainly focus on compatibility of legacy json (only for internal testers, we havenot published it).
I think json file specs are scattered everywhere in this file. That is a mess (blame GPT for it's low context window for 20USD user)

# _load_ion_module
## line  223 - 253
@Henry: Called in try_merge_selected_sample and run_pseudolabeling for now.
It aims to find out possible important fragments(peaks) in spectral data
I think it is bind to data preprocessing workflow (Under button `Prepare Dataset`)
Do we really need to import here, not when user clicks the Prepare Dataset button? (Second-level operation)

# Changelogs
## line 256 - 336
@Henry: need to be moved somewhere in future

# tkinter-safe init
## line 338 - 345
@Henry: remind me of CLI in future

# raw/mzml file conversion init
## line 347 - 348
@Henry: I should think about moving file conversion into another sub-window as data preparation and ML analysis.
The init here makes the code not elegant, like a summer project

# _to_posix and _from_posix
## line 352 - 353 (removed, no longer exist)
@Henry: the questioned comment is correct. Claude Code confirmed pathcanon does that work well
Replaced _from_posix with str(pathcanon.to_native_path(p))

# pure definition: FILETYPE_TO_KEY
## line 350 - 360
@Henry: Called in later drag and right click events (around line 5500)
Should this be applied by default, or just specifically to drag and right click?

# integrity check (placeholder)
## line 361 - 374
@Henry: Placeholder. Expected to work in future when the files/projects are migrated to different devices/paths
If Claude Code has the ability to do this, draft the spec to me and I will see if we can make it work soon


# App logger 
## line 376 - 406
@Henry: Shows log in main GUI loop and save button for dumping all logs presented there.
simple debug print won't trigger this, if we need, use logger.log (defined in line 405)

# Class:PseudolabelingSetupWindow -> CGASetupWindow
## line 408 - 1071
**Include functions inside, quick check is performed**
@Henry: replace with CGASetupWindow. Lots of hard-coded definition, needs carefully maintanience
(Probably need a doc to point the place to modify)

# Class:MetadataEditorWindow
## line 1074 - 1600
**Include functions inside, quick check is performed**
@GPT: add a bug tracking print to identify converted csv not finalized issue on internal testers [BUG detection]
@Claude Code: remove cwd guess by abs path (root cause of error emerges)


# write_to_gui
## line 1602 - 1609
@Henry: asked Claude Code of why it is here, and noticed tkinter is chained here to AppLogger.
This means original AppLogger supports non-tkinter calls. This is intended during design phase,
probably when I got error testing CLI functions and asked GPT of why AppLogger throws tkinter errors.

# `get_version_info_for `
## line 1611 - 1623 
@Henry: Kept.
@Henry: Exactly analysis from Claude Code is correct. Maybe change the role of this function in future
Thanks. I will mark get_version_info_for as never called and consider other  
  uses of it ()


# Data Loading helpers
## line 1625 - 1713
@Henry: confirmed that robust_read_csv is a well-organized csv/tsv handler in GlycoMSP, and should be 
widely called except for specific (limited) usage. Replaced 2 (_read_any_table and extract_ion_intensities)
`_robust_read_csv` is changed to `robust_read_csv` since it is generic function

# Main GUI - About Window
## line 1715 - 1752
@Henry: About window. Fixed the path incorrected loaded (not speficied in configparser, causing cwd path)issue

# `validate_file_path`
## line 1754 - 1763
# Kept as it-is
@Henry: That was called when selecting Raw file/mzML in main window. The analysis was incorrect.
tk.Button(button_frame, text="Select Raw File", command=lambda: select_file("raw")).grid(row=0, column=0, padx=5)
select_file calls validate_file_path inside it
(I think the reverse-engineering would reveal this conflict, maybe)

# Main GUI-related action
## line 1764 - 1779
@Henry: update_display, clear_files and reset_main_status are located here. No more?

# build_final_path (duplicates) and finalize_converted_files
## removed
@Henry: seems weird and found better work in generate(self) from MetadataEditorWindow. Confirmed with Claude Code.

# select_files_generic
# removed
@Henry: asked Claude Code and noticed the real calling function (identical) in open_prepare_dataset_window

# `_norm_for_display`
# Removed as suggested.
@Henry: I think that's a surgical fix for tree not displaying properly.

# `create_tkinter_root (line 1738)`
# Removed
@Henry: Why it is still there? Is that the bug source of incorrect Window priority (level)?
**Blocks removed**
# Function to create and hide Tkinter root window
# If future version supports GUI fully available, please reconstruct this part (not hiding main window, but think about what should be there)
def create_tkinter_root():
    root = tk.Tk()
    root.withdraw() # hide main window
    return root
**removed **

# select_file for spectral file selection in Main GUI
## line 1781 - 1795
@Henry: updated to raw/mzml only


# same_drive and promote_temp_to_final
## line 1830 - 1857
@Henry: Moved next to select_file since these 2 are function modules for ms2/ms3 csv rename & finalization


# on_batch_conversion_complete
## removed
@Henry: asked Claude Code and confirmed the behavior has been integrated in MetadataEditorWindow

# launch_metadata_for_all
## removed
@Henry: called from on_batch_conversion_complete so remove as well
@Claude Code: ⏺ Yes — MetadataEditorWindow handles the entire batch internally. Here's the    
  loop:                                                                    
  1. Constructor (line 1090): receives raw_file_list (the full batch), sets     
  current_index = 0, immediately calls self.load_file(self.raw_file_list[0])    
  (line 1125).                                                                  
  2. load_file → starts background_conversion in a thread → on success schedules
   on_conversion_complete → enables metadata fields for user input.             
  3. User clicks "Generate Dataset" → generate() runs → saves metadata JSON,
  renames temp CSVs, fires self.callback(raw_path, metadata, savename).         
  4. Batch advance (lines 1585-1593): increments current_index. If more files
  remain, calls self.load_file(next_file) — the window reuses itself for the    
  next file. If last file, shows completion dialog, destroys window, calls
  self.on_finish().                                                             
                                                            
  So the batch loop is: __init__ → load_file → generate → load_file → generate →
   ... → destroy. The window processes files one at a time, reusing itself
  across the batch. launch_metadata_for_all would do the exact same thing —     
  create a MetadataEditorWindow with a file list and let the class handle the
  loop. It's completely redundant with what launch_metadata_batch already does
  at lines 1899 and 1943.

  Safe to remove launch_metadata_for_all and on_batch_conversion_complete.  

# Metadata-related supporter block
# on_metadata_ready
## line 1860 - 1999
@Henry: asked Claude Code to improve the log output.
@Claude Code: metadata is a full dict — not useful in a one-liner log. Leave it unused but  
  in the signature. It's there for future use (e.g., auto-assigning to experiment_projects).     

# _score_metadata_candidate and _choose_metadata_dialog
# removed
@Henry: asked Claude Code and confirmed it was fuzzy logic metadata resolution. 
Now we use a strict file validation as our plan, so this logic is not longer applicable.

# _fallback_simple_ion_scoring
## line 2001 - 2032
@Henry: a safe net confirmed by Claude Code when marker ion hit list hits exception.

# CGA_glycotope_score
## removed
@Henry: triple-quoted string with copy-paste of _fallback_simple_ion_scoring, probably placeholder for CGA score B.
No need to keep since we have real CGA-B (although it is workbook-based)

# UID related helpers: _infer_uid_series _ensure_uid_and_origin and _outer_union_concat
## line 2037 - 2093
@Henry: asked Claude Code and it said The curated library export
(msp_curated_library_db.export_trainable_csv) will also need to combine datasets with UID tracking.
It was prepared for later combine_trainable_datasets_ui to keep combined dataset trackable.
Notice that the breakpoint of trainable file isolated from current method json. 
Until we include ML method recorded and revealed in prepare dataset window, these are orphans in current version.

# combine_trainable_datasets_ui
## line 2095 - 2174
@Henry: Combine Trainable Datasets called in ML analysis window


# open_prepare_dataset_window
## Window of Prepare dataset including GUI and widgets inside, REALLY LARGE
## line 2176 - 6000
@Henry: Will decompose modules inside. Too big to summarize in one line.

# path resolver: _norm(p) and _rel
## removed
@Henry: Since it was never called, asked for investigation and Claude Code said that,
this is the work of converting file to relative path. If we need portable one-click file conversion, add then.

# path resolver: _abs
## line 2210 - 2212
@Claude Code:Actively used to resolve POSIX paths from experiment JSON back to absolute native paths.  

# _get_selected_exp_sample
## removed
@Henry: confirmed with Claude Code that more flexible _get_selected_context is working version to accept CGA score B
and any future file entry expansion

# _current_exp_title
## line 2212 - 2228
@Claude Code:For passing clean Experiment name. Called in _get_ml_context (not used) and save_current_experiment_method

# _infer_method_kinds
## line 2232 - 2241
@Henry: Legacy recognition of MAS/CGA method file, for later exp json save.
_detect_legacy_family doesn't require csv (not for export/saving exp json) so can't be replaced


# _ensure_method_stub (shadowed)
## removed 
@Henry: asked Claude Code and the later one is safer (handles None) and it is called from around line 6000~
so the early definition around 2600 is fine even we set the working one around line 2900

# _ensure_methods_before_saving_exp_v1
## line 2244 - 2303
@Henry: force method file appears before saving exp json. Produce one method entry per family

# exp.json save/load helpers block
## line 2305 - 2610
Contains import_experiment_json, export_experiment_json and removed functions below:
##_infer_family
## removed
@Claude Code: Recommendation: Replace with _detect_legacy_family(files) or "UNKNOWN" and delete    
  _infer_family. One fewer near-duplicate.
## import_experiment_json_old
## removed
@Claude Code:The new version (line 2648) is a strict superset


# integrity check placeholder
## line 2612 - 2628 
@Henry: future work in v1.1 or v1.2

# _get_selected_context
## line 2631 - 2700
@Henry: confirmed with Claude Code that this is clean and working one 

# _ensure_method_stub (active one)
## line 2702 - 2726
@Henry: ensure method node exists, if no, add an unknown for it

# _detect_legacy_family
## line 2728 - 2741
@Henry: the none-asking csv version of _infer_method_kinds. Legacy method load use.



# _get_selected_context_old
## removed
@Henry: use the new one, but note the changes to new one
@Claude Code: The active version is a cleanup of the old one — same        
  algorithm, same tree-climbing approach, same method resolution logic. Changes: (1)
  removed redundant exp_node null guard, (2) consolidated files lookup to one place    
  instead of two, (3) combined exp name parsing into one line. No behavioral
  difference. The old version was the draft, the new one is the polished copy.
  
# _register_or_update_method_ref
## line 2743 - 2770
@Henry: main updater for method (called when method file needs to be updated).

# _export_method_v1
## line 2772 - 2808
@Henry: main method export module.
Removed v1 in filename (we have it in schema) to avoid file selection not recognising by extension (.method.json)

# definition area of parameters in sub-windows from Prepare Dataset Window
## line 2810 - 2830
@Henry: asked Claude Code if I can move it, if we do so, need to put them before getting called.
Current structured format is also good.

# Treeview: data explorer panel
## line 2832 - 2835 
@Claude Code: #0 — Treeview's built-in first column. Every Treeview has a hidden column called #0
   that displays the tree hierarchy (expand/collapse arrows + text). You can add       
  additional named columns via columns=("col1", "col2"), but this tree only uses #0 —  
  it's a tree-only display with no extra data columns.
 
# Tree logic block

# on_tree_select(event)
## line 2838 - 2856
@Claude Code: Updates the Link and Merge button states based on what the user clicked in the tree.
event: Keep it. It's required by tkinter's .bind() API. Line 3182 binds this   
  function via tree.bind("<<TreeviewSelect>>", on_tree_select) — tkinter always passes
  an event object to bound callbacks. If you remove the parameter, you get a TypeError 
  at runtime. The function doesn't need to read anything from event (the selection is
  accessed via tree.selection() instead), but the parameter must stay in the signature.
@Henry: Check validation when select a tree node (select sample within a method should be clickable in this version).


# tree.bind("<<TreeviewSelect>>", on_tree_select)
## line 2858 - 2859
@Henry: Tree-view selection event listener(handler), calls on_tree_select when elements in panel get clicked

# try_link_selected_sample
## line 2861 - 2881
@Henry: Link & Validation on MAS method sample. Considering CGA compatibility
removed leftover line

# try_link_selected_sample_old
## removed
@Henry: asked Claude Code and confirm older version that only works, 
if the user clicks exactly on the sample node.

# refresh_tree
## line 2883 - 3007 
@Henry: Feeling the selection drifting bug coming here.
@Claude Code: 
  - Block A (lines 3249-3257) — Sample-level anchors: Shows CSV, Excel, and Metadata   
  directly under the sample node. These are "identity" files shared across all methods
  — they don't belong under any specific method.                                       
  - Block B (lines 3259-3272) — Build method list: Reads files["_methods"] if it
  exists. If not, falls back to constructing a single-entry list from files["json"]    
  (legacy compat). Produces a methods list for block C to iterate.
  - Block C (lines 3274-3331) — Insert method nodes: For each method, creates a Method 
  node with family tag. Under it: the method JSON path, then CGA artifacts (insilico,  
  ionlist, Score B, pseudolabel) under idx == 1 only (Phase 1.5 limitation — artifacts
  aren't stored per-method yet). For MAS methods, reads the method JSON file directly  
  to show Converted CSV / Annotation Excel / Metadata from its inputs dict.
@Claude Code:
  The symptom you describe — "clicking sample A, then selecting sample B, and the      
  previous entry shifts" — matches this pattern. Every refresh_tree() call wipes and
  rebuilds, and if something triggers refresh_tree between your click and your next    
  action (e.g., a background callback, a method save), the selection jumps.
  A fix would be: Save the selection before rebuild, restore it after
@Henry: Fix applied. In future we may need to add sample node edit/rename/deletion(now the delete is buggy)
May need redesign then (v1.1 or v1.2+)
@Claude Code: the item_id

# show_experiment_summary
## removed
@Henry: reason - we have icon showing already. Not worth playing math here.

# clean_sample_name
## line 3010 - 3028
@Henry: Avoid status symbol contaminate the tree contexts
@Claude Code: No risk of undefined error. All callers (lines 2881, 2913, 3055, 3204, etc.) are     
  inside nested functions that only execute at runtime (button clicks, tree events,
  callbacks). None of them run at function definition time. By the time any of them    
  fires, Python has already executed past line 3372 where clean_sample_name is defined.
  
# add_sample():
## line 3030 - 3052
@Henry: Add sample (customized empty node) to selected experiment
Caveats: only can add new empty sample node when experiment already exists.
Caveats: can't add same sample name (future plan: allow rename/auto-numbering) after multi-method on single sample supported


# minor definition fix
## 20260415 (3 places)
experiment_projects[exp_name]["samples"][sample_name] = {"csv": None, "excel": None, "json": None}
experiment_projects[exp_name] (or [exp_title])["samples"][sample_name] = {"csv": None, "excel": None, "metadata": None, "json": None}


# clean_unassigned_samples
## line 3054 - 3071
@Henry: Temp stable solution to remove empty sample or experiment with no sample inside.

# store_metadata_back_to_sample
## line 3073 - 3113
@Henry: chained with open_metadata_editor_for_sample that opens MetadataEditorWindow to create post-metadata
 when user clicks Link Sample button on a sample with both CSV and EXCEL (MAS workflow)
Updates: json_path -> meta_path (readability), remove legacy write_method_file line,
and add a fallback that allow link and validation continues after metadata is created.

# open_metadata_editor_for_sample
## line 3115 - 3147
@Claude Code: before your __init__ patch, this was the workaround for the same problem:    
  __init__ called load_file too early, so this site used after(10) to delay it. Now    
  that __init__ handles the ordering correctly, this workaround is unnecessary. Safe to
   comment out.

# Negative Sampling Options block   
## line 3149 - 3471 

# current_selected_files
## line 3150 - 3158 
@Claude Code: resolves the tree selection to (exp_name,            
  sample_name) via _get_selected_context(), then returns the files dict for that sample
   from experiment_projects.
   
# current_selected_files_old
## removed
@Claude Code:  The new version (line 3506) just calls _get_selected_context() and looks up the      
  result in experiment_projects. Three lines vs eighteen lines, same result. The
  refactoring happened when _get_selected_context was introduced as the centralized    
  selection resolver. The old version also has a bug... (skipped in this report)
  
# _fetch_ion_count_for_dialog
## line 3161 - 3184
@Henry: Negative Sampling Excel picker

# _update_gate_hint
## line 3186 - 3197
@Henry: Live update the tk var value in Negative Sampling Options when the ion(feature) counts are known

# open_negative_options_dialog
## line 3199 - 3471
@Henry: Prepare Negative Sampling Options Panel that passes/update the negative sampling accross MAS (merge MAS) and CGA (CGA to trainable) 

# Negative Sampling Options block ends


# Ion Mining block
# open_ion_suggest_dialog -> open_ion_mining_dialog
## line 3273 - 3302
@Henry: name changed to meet manuscript description

# Link and Validate sample before merge in MAS workflow

# link_and_validate_sample
## line 3304 - 3395
@Henry: remove old commented lines, refine the logger message, and support fragmentation_mass (mass) in score B workbook-based

# _resolve_method_json
## removed
@Henry: I think it is not called and the only site calling it might be a mistake (see try_merge_selected_sample)
@Claude Code: The filesystem search strategy it offered (CSV folder, Excel folder, CWD) was a      
  workaround for stale/relative files["json"] paths. But _get_selected_context already 
  resolves method paths via the _methods list, and metadata paths come from            
  files["metadata"] directly. No remaining code needs the fallback search.

# _get_ml_context
## removed 
@Henry: either exp name (current_selected_files) or method file can be retrieved from _get_selected_context


# try_merge_selected_sample
## line 3397 - 3525 
@Henry: Identified a previous bug of referencing method json and asking derivatization field wrongly.
Proposed code fix can be found in Claude Code session logs. Aggressive test only confirms _resolve_method_json
and the method json is correctly referenced and the asking derivatization returns "Others"
Added extra tracker:
            # 20260416 add trainable_csv tracking (in future we may need subnode, or separated ML method under same MAS/CGA method)
            files["trainable_csv"] = outpath 
In CGA -> trainable I can't find the entry, need discussion later


# MAS to trainable merge workflow block ends

# CGA analysis result (previously called Pseudolabel) data to trainable csv functionality block

# _canonicalize_scan_column
## removed
@Henry: feel like _ensure_scan does better.
@Claude Code: _ensure_scan is a superset of _canonicalize_scan_column with fallback and expection handling.


## _coerce_comp_to_tuple
# removed
@Henry: calls _TUPLE_LIKE_RE (re.compile(r"^\(\s*\d+\s*,\s*\d+\s*,\s*\d+\s*,\s*\d+\s*,\s*\d+\s*,\s*\d+\s*\)$"))
to produce pure numeric composition tuple such as (6,5,3,2,1,0)
Not called in current version ,probably the bridging of tuple input has been fixed previously
@Claude Code: If you ever need string→tuple parsing, it's a trivial ast.literal_eval one-liner.

# SCAN_CANDIDATES = ("MS2scan_no", "MS2Scan_no", "ScanNum", "scan", "Scan", "unique_ID")
## removed and replaced duplicated definition

# _NONMASS
## removed after clarifying _infer_ion_masses_from_wide_df (only parse trainable csv features for conversion)

# _infer_ion_masses_from_wide_df
## line 3539 - 3551
@Henry: asked Claude Code and trace the origin, should be something reusing trainable csv features
The word choice in that panel was a disaster, either change the description or name.

# _ensure_scan
## line 3553 - 3579
@Henry: return out = df.copy() with a guarantee of canonical string column 'MS2scan_no' exists in df, 



# _comp_tuple_to_label
## line 3589 - 3626
@Henry: asked several questions, added non-interget ValueError warns. 
@Claude Code:
⏺ Question 1: int(float(p)) — The float() first handles cases where the string contains
   decimal values like "2.0" from CSV serialization. int("2.0") raises ValueError, but
  int(float("2.0")) works → 2. The risk is truncation of actual fractional values      
  (e.g., "1.5" → 1), but monosaccharide counts are always integers, so this is safe for
   the domain.                                                                         
                                                                                       
  Question 2: You're right — there's no input validation of tuple order. The function  
  assumes the input tuple is always (Hex, HexNAc, NeuAc, NeuGc, KDN, Fuc) (line 4093
  comment). The mapping at line 4097 destructures in that exact order. If upstream ever
   produces a different order (e.g., Fuc-first like (Fuc, Hex, HexNAc, NeuAc, NeuGc, 
  KDN)), the output labels would silently assign wrong counts to wrong sugars. This
  order is set by compnewv4 — as long as the composition enumerator doesn't change its
  tuple order, it works. But it's a fragile implicit contract with no runtime check.

  Question 3: if len(x) < 6: x = list(x) + [0]*(6-len(x)) — Zero-pads short tuples to  
  length 6. If someone provides (5, 4, 2) (only 3 elements), it becomes (5, 4, 2, 0, 0,
   0) so the destructuring at line 4097 doesn't crash with "not enough values to       
  unpack." It assumes missing trailing positions are zero — which means missing values
  are treated as NeuGc=0, KDN=0, Fuc=0. This is reasonable for N-glycan compositions
  where KDN and NeuGc are often absent.

# _count_hits_to_ionlist
## line 3628 - 3666
@Henry: function that provides minimal hit count for negative sampling: decide lower boundary of glycan vs Non-glycans
@Claude Code:  _to_array — Converts peaklist/peakintensity from whatever format they come in     
  (list, tuple, JSON string "[100.5, 200.3]", Python literal string "(100.5, 200.3)",  
  or comma/semicolon-separated string "100.5;200.3") into a numpy float array. This    
  handles the inconsistency where CSVs store peaks as strings. 


# build_trainable_from_pseudolabels block (CGA/PL -> trainable) -> build_trainable_from_CGA
@Henry: change name to CGA                       
# single line log = (logger.log if logger else print)
## line 3696
@Claude Code:If a logger object was passed in (the AppLogger instance), log becomes logger.log. 
If logger is None (e.g., called from CLI or tests), log falls back to print.
@Henry: In future CLI support we may need this implementation globally

# single line aliases = {c.lower(): c for c in pl.columns} 
@Claude Code: It builds a dict mapping lowercased column names to their         
  original-case names. For example if pl has columns ["MS2scan_no", "PeakList",        
  "Composition"], aliases becomes:                                                     
  {"ms2scan_no": "MS2scan_no", "peaklist": "PeakList", "composition": "Composition"}                                                                                      
  Then has_lower("peaklist") checks "peaklist" in aliases → True, and                  
  aliases["peaklist"] gives back "PeakList" — the actual column name to use with the   
  DataFrame. It's a case-insensitive column lookup mechanism, because different        
  upstream tools may write the same column with different casing.

# build_trainable_from_pseudolabels -> build_trainable_from_CGA
## line 3667 - 4260
@Henry: Found critical bugs, fixed :reuse branch is broken (halfway refactored copypaste code from rebuild)
Change rebuild name to extract and update description. Fixed and started sharing workflow of 2 method after ion_masses assignment.
In stage 5 (add negatives) a possible defense point 
(specifically ask cowork to extract this with the stage summary,
 and replace the block with a short summary telling Claude Code the tech details and manuscript concern moved to another file.
 You create one markdown file for saving that):

❯ Finished reading stage 5 block. I see a possible logical issue after our discussion  
  on stage 4 last night. In stage 5 we don't do a topn selection, and for same         
  MS2scan_no assigned as different compositions AND having low ion score that falling  
  to negatives will be kept. I took a look at PL/CGA tsv and converted trainable csv   
  and confirmed this suspect. We have plenty of negative labels for downstream machine 
   learning process, so this really doesn't matter that much. What do you think?       
  Should we fix this bug by introducing similar topn selection in stage 4, or just     
  note it for now?                                                                     
                                                                                       
⏺ Note it for now. The duplicates slightly inflate the negative count per scan, but    
  since negatives are already capped by max_ratio (line ~4712: keep =                
  min(len(neg_rows), int(max_ratio * max(1, n_pos)))), the excess gets trimmed anyway. 
  And ML-wise, having a few extra negatives from the same scan with different failed
  compositions doesn't hurt — they all represent "this scan didn't match glycan ions   
  well," which is the correct signal.                       

  Fixing it would add complexity (dedup by scan before or after the hit count loop) for
   minimal ML benefit. If you ever see class imbalance issues traced to this, revisit
  it then.                                                                             
                                                            
❯ Got it. Then if reviewer ask me why duplicated MS2 scan with negative label presents 
   in the trainable dataset, I know the reason and can defend by cap and nearly        
  non-impact outcome. FYI since those are labeled as Non-glycans (composition          
  assignment rejected) now, it becomes simple duplicated entries in trainable data.    
                                                                                       
⏺ Exactly. And if a reviewer presses further, you can add: the duplicates are from     
  different failed composition assignments for the same scan — they share the same peak
   data but were independently rejected by the ion score threshold. After feature      
  extraction, they produce identical feature vectors with the same "Non-glycan" label,
  so they're effectively duplicate training rows. Duplicates in the majority class get 
  capped by max_ratio balancing anyway, making the ML impact negligible. Clean defense.







⏺ Here's the stage summary:                                                            
                                   
  Stage: 1) Load & prep                                                                
  Lines: 4177-4205                                                                     
  Operation: Load pseudolabel TSV, clean columns, find scan column, resolve aliases    
  Key inputs → outputs: pseudo_path → pl DataFrame, pl_scan_col                        
  ────────────────────────────────────────                                             
  Stage: 1b) MS1 attach                                                                
  Lines: 4208-4244                                                                     
  Operation: If mass feature or precursor gate needed, merge protonatedmass from       
    converted CSV into pl                                                              
  Key inputs → outputs: converted_csv_path → pl gains protonatedmass column         
  ────────────────────────────────────────                                             
  Stage: 1c) Precursor gate                                                         
  Lines: 4230-4244                                                                  
  Operation: Filter rows by PPM distance between observed and theoretical mass         
  Key inputs → outputs: pl filtered by precursor_gate_ppm                           
  ────────────────────────────────────────                                             
  Stage: 1d) Composition normalize                                                  
  Lines: 4252-4272                                                                     
  Operation: Ensure composition column exists (from comp_str, comp_tuple, pseudo    
    compositions)                                                                      
  Key inputs → outputs: pl gains canonical composition column                       
  ────────────────────────────────────────
  Stage: 2) Thresholding / selection                                                   
  Lines: 4276-4414
  Operation: Filter by min ion score + max PPM, take top-N per scan. If Score B columns
                                                            
    exist: prefer score_b_selected rows, fall back to Score A for uncovered scans.
    Normalize Structure column.
  Key inputs → outputs: pl → sel (selected positive rows)
  ────────────────────────────────────────
  Stage: 3) Ion list                                                                   
  Lines: 4419-4423
  Operation: Load ion list file                                                        
  Key inputs → outputs: ion_file_path → ion_df              
  ────────────────────────────────────────
  Stage: 4a) Features — REBUILD                                                        
  Lines: 4431-4571
  Operation: Extract features from raw peaks via build_features_from_peaks_log10_plus1,
                                                            
    merge Structure labels, attach MS1 mass if requested, handle NaN Structure cleanup
  Key inputs → outputs: sel + peaks → pos_feat with features + labels
  ────────────────────────────────────────
  Stage: 4b) Features — REUSE                                                          
  Lines: 4573-4673
  Operation: Load pre-built wide CSV, infer ion masses from headers, merge labels      
  Key inputs → outputs: wide_feat_csv → pos_feat with features + labels
  ────────────────────────────────────────
  Stage: 5) Negatives                                                                  
  Lines: 4675-4731
  Operation: Scan PL for rows NOT in positive set with < min_hits ion matches. Cap by  
    ratio. Build features for negatives. Label as "Non-glycan"
  Key inputs → outputs: pl scans not in pos_ids_set → neg_feat
  ────────────────────────────────────────
  Stage: 6) Finalize + save                                                            
  Lines: 4732-4811
  Operation: Concatenate pos + neg, ensure scan column, drop helper columns, order     
    columns [MS2scan_no, Structure, ...features], save CSV, build summary dict
  Key inputs → outputs: pos_feat + neg_feat → final_df → output_path


# ion_massed investigation
@Claude Code:  ⏺ Here's the summary:                                                                   
  Line: 4431                                                                           
  Stage: 4) Features — REBUILD mode                                                    
  Source: ion_df["mass"] from ion list file                                            
  What happens: Reads masses from the loaded ion list DataFrame. Clean, explicit       
  source.                                                                              
  ────────────────────────────────────────                                             
  Line: 4435                                                                        
  Stage: ↳ assignment                                                                  
  Source: ion_df["mass"].astype(float).tolist()                                        
  What happens: ion_masses = list of floats from ion list
  ────────────────────────────────────────                                             
  Line: 4447                                                                           
  Stage: ↳ usage
  Source: passed to _count_hits_to_ionlist(...)                                        
  What happens: Counts ion hits per scan for Score A        
  ────────────────────────────────────────
  Line: 4462                                                                           
  Stage: ↳ debug
  Source: len(ion_masses)                                                              
  What happens: Logging                                     
  ────────────────────────────────────────
  Line: 4589                                                                           
  Stage: 4) Features — REUSE mode
  Source: _infer_ion_masses_from_wide_df(feat)                                         
  What happens: The risky one — guesses masses from column headers of pre-built wide
    feature CSV
  ────────────────────────────────────────
  Line: 4690                                                                           
  Stage: 5) Negatives
  Source: pd.DataFrame({"mass": ion_masses})                                           
  What happens: Reconstructs ion_df from ion_masses if original ion_df is empty
  ────────────────────────────────────────
  Line: 4724                                                                           
  Stage: ↳ usage
  Source: build_features_from_peaks_log10_plus1(neg_df, ion_masses, ...)               
  What happens: Builds features for negative samples        
  ────────────────────────────────────────
  Line: 4807                                                                           
  Stage: 6) Finalize
  Source: len(ion_masses or [])                                                        
  What happens: Logged in output summary                    
  ────────────────────────────────────────
  Line: 8679/8692                                                                      
  Stage: try_merge_selected_sample (separate function)
  Source: passed through from merge flow                                               
  What happens: Feature building for Route B merge          



# open_ion_suggestions_viewer
## line 4263 - 4340
@Henry: just a sub window of viewing csv from ion mining suggestion. Safe and don't touch until user wants improvements

# Run CGA analysis block (previously called Pseudolabeling)

# append_runlog
## line 4343 - 4364 
@Henry:Create jsonl as run log when a CGA analysis is completed (generates scored A + B(optionally) tsv)

# run_pseudolabeling -> run CGAanalysis (in next version, not applied yet)
## line 4366 - 4720



# comapare lib vs libn
@Claude Code
⏺ lib — the raw in-silico CSV as loaded from disk. Contains columns like Hex, HexNAc,  
  NeuAc, NeuGc, KDN, Fuc, Mass, and possibly HexA, SO3, PO3H if the library was        
  generated with modifiers.                                                            
   
  libn — the normalized version after marker.normalize_insilico(). The normalizer      
  standardizes the library: creates canonical comp_tuple (6-element tuple) and comp_str
   (string representation) from the base composition columns, and ensures Mass is      
  sorted and consistent. But the normalizer only knows about the 6 base columns — it
  doesn't carry through the modifier columns (HexA, SO3, PO3H).

  libn[_c] = lib[_c].values — This re-attaches the modifier columns from the original  
  lib onto the normalized libn. .values extracts a numpy array, bypassing pandas index
  alignment. This works because normalize_insilico preserves row order and count       
  (one-to-one mapping) — row 0 in libn corresponds to row 0 in lib. The .values is
  important: without it, if lib and libn have different pandas indices (e.g., after
  normalize_insilico reset the index), a direct libn[_c] = lib[_c] would align by index
   and potentially produce NaN mismatches. Using .values forces positional assignment.
   
# run_score_b_enrichment
## line 4722 - 4767
@Henry: the funtion called by start CGA analysis (automatically)
when score_b_workbook_path exists, score A and score B will be both calculated.
Wait, where is the rescue score B from a CGA-analysis result file? We don't have it?


# assign_file
## line 4769 - 4776
@Claude Code:takes a file type key, a file path, an    
  experiment name, and a sample name, ensures the experiment and sample containers
  exist in experiment_projects, sets the file path, and refreshes the tree
  
# extract_rawname_from_metadata
## line 4780 - 4792
@Claude Code: Silently returns None for any failure — malformed JSON, wrong type, missing field — 
so the caller falls back to "Unassigned".

# extract_title_from_metadata
## line 4794 - 4804
@Henry: moved and grouped with extract_rawname_from_metadata as metadata extractor block

# handle_csv_selection
## line 4809 - 4812

# handle_excel_selection
## line 4814 - 4817

# handle_json_selection
## removed
@Henry: replaced with next version: handle_metadata_selection

# handle_metadata_selection
## line 4819 - 4836
@Henry: removed unused sample_name

# select_files_generic
## line 4840 - 4868
@Henry: expand filetype dict and add crash-preventing try-exception condition

# update_status_display
## line 4871 - 4880
@Claude Code:It's not useless by design — it just has a narrow scope (only works for experiments  
  whose status labels were registered). Keep it, but it's low-value in current UX since
   users mostly work under named experiments, not "Unassigned." When you redesign the  
  status bar in a future pass, this function would be the place to expand.

# TREEVIEW OPERATION block

# on_drag_start
## line 4882 - 4905
@Henry: asked Claude Code for debug print suggestion, now print won't show iid (for tkinter reference, not fixed)
@Claude Code: debug print line is located at wrong position (suggested fix applied)

# on_drag_release
## line 4907 - 4931
@Henry: removed GPT symbol enthusiast condition check

# on_right_click
## line 4934 - 4985
@Henry: fixed filetype and filename extra strip (corruption) issue and link sample if statement condition bug

# move_file
## line 4987 - 4992
@Claude Code: Three steps:
  1. Grab the file path from the source sample's files dict using the file type key    
  (e.g., files["csv"] → "/path/to/data.csv")                
  2. Verify the basename matches what the user sees in the tree (filename from the tree
   label text). This prevents stale moves — if the tree shows "data.csv" but the dict  
  was updated since the menu was built, the basenames won't match and the move silently
   aborts.                                                                             
  3. Move it: Null out the source (= None), then assign_file puts the same path into
  the destination experiment/sample. The actual file on disk doesn't move — only the   
  reference in experiment_projects changes. assign_file also calls refresh_tree() at
  the end.                                                                             
  It's not magic — it's just a reference reassignment. The file stays on disk where it 
  is; only the in-memory dict bookkeeping changes. Think of it as "unlinking from
  sample A, linking to sample B."
 
 
 # discussion about auto-update method file
@Claude Code: For auto-update in the future, save_method_v1_for_sample with auto=True already      
  serves that role — it resolves the output path automatically without dialogs. No need
   for a separate function. 
   
# _save_method_v1_now
## removed
@Claude Code: Uncalled and the role is covered by other functions

# save_method_v1_for_sample
## line 4998 - 5061
## line moved to near _export_method_v1
@Henry: ask Claude Code to explain multiple times.
@Claude Code:
⏺ We already covered save_method_v1_for_sample in detail earlier — it builds a method  
  v1 dict from the current tree state (via build_method_v1_from_tree), writes it to    
  disk as JSON, and registers it in files["_methods"] via                              
  _register_or_update_method_ref. It's the core method-save operation called from
  multiple places. 
  About moving it: yes, move it near _export_method_v1. Both are method JSON writing   
  functions — _export_method_v1 is the UI-triggered wrapper (button click),
  save_method_v1_for_sample is the programmatic backend. They belong together          
  logically. Right now save_method_v1_for_sample is stranded between move_file and the
  tree UI wiring code, which makes no sense organizationally.

# write_method_file
## line 5063 - 5136
## line moved to near _export_method_v1
@Henry: Legacy writing fallback only (final failsafe to prevent method file is not created)
@Claude Code: If your test confirmed no regression, delete it entirely rather than just moving it.


# safe_get_field
## deleted
@Claude Code: Every place in v14 that needs this pattern already uses .get() directly. Delete it.

# _resolve_path
## line 5139 - 5155
@Claude Code:
  1. Empty/None/non-string → returns None (lines 5769-5775)                            
  2. Already absolute (POSIX /path/... or Windows C:\path\...) → returns as-is (lines
  5778-5783)                                                                           
  3. Relative path → joins with base and normalizes:        
  os.path.normpath(os.path.join(base, p)) (line 5785)                                  
                                                            
  The Windows drive check at line 5782 is important — on macOS/Linux,                  
  os.path.isabs("C:\\Users\\...") returns False because it doesn't recognize Windows
  drive letters. The regex catch ensures Windows absolute paths aren't accidentally    
  treated as relative and joined with the base directory.   

  This is the core path resolver for normalize_method_json — every file path stored in 
  a method JSON goes through this to convert from however-it-was-stored to an absolute
  path usable at runtime. It's essentially what _abs does (the function we kept        
  earlier), but with the extra Windows-on-macOS safety net. 
  
# normalize_method_json
## line 5158 - 5339
@Henry: replace pseudolabels_tsv -> CGAresult_tsv

# build_method_v1_from_tree
## line 5341 - 5463
@Henry: replace pseudolabels_tsv -> CGAresult_tsv
 
# load_method_file
## line 5465 - 5520
@Henry: added some log outputs

 
# change_sample_method_folder
## line 5523 - 5532
@Henry: asked Claude Code and get fix of folder path is always None issue (global reference None on the top)
@Claude Code:
⏺ Here's the full flow:                                                                
                                                                                       
  sample_method_folder is a closure variable initialized to None at the top of         
  open_prepare_dataset_window. It stores the directory where per-sample method JSONs
  get auto-saved. Functions like save_method_v1_for_sample (line ~5616) read it:       
                                                            
  base_dir = sample_method_folder or (os.path.dirname(sample["csv"]) if
  sample.get("csv") else os.getcwd())                                                  
   
  If sample_method_folder is set, method JSONs go there. If not, they fall back to the 
  CSV folder.                                               
                                                                                       
  change_sample_method_folder lets the user pick a directory via dialog and stores it  
  in that variable. But with the current global declaration, it writes to module scope
  — the closure variable stays None. That's the scoping bug.                           
                                                            
  update_status_display reads sample_method_folder directly (line ~5446):              
   
  sample_label.set(f"Sample method will be saved at: {sample_method_folder if          
  sample_method_folder else 'None'}")                                                  
   
  With the global bug, this always shows "None" because the closure variable was never 
  updated. After fixing to nonlocal, the user picks a folder → sample_method_folder
  updates in the closure → update_status_display reads the updated value → status bar  
  shows the actual path.                                    

  So the fix connects the full chain: button click → folder stored → status bar        
  reflects it → method saves go to the right place.

# change_sample_method_folder
## line 5534 - 5539
@Henry: Claude Code fixed the global bug making it always None (on the top of global definition)

# load_experiment_method_file
## line 5541 - 5554
@Henry: Claude Code fixed a Error (Why I didn't get stopped from loading exp json? Is that silent bug?)

# save_current_experiment_method
## line 5556 - 5578
@Henry: Claude Code fixed the behavior of automatically point to current experiment json saving folder

# bind_right_click
## line 5580 - 5585
@Henry: Claude Code suggested a fix on MacOS right-click not effective issue

# get_selected_csv_path
## removed
@Claude Code: This was an early prototype that misunderstood the data structure.
Delete it. current_selected_files() is the proper way to get file paths from the selection.

# launch_pseudo_labeling
## line 5587 - 5866
@Claude Code: 
⏺ Your understanding is correct — launch_pseudo_labeling is the orchestrator that:     
  1. Resolves the selected sample's context (files, metadata)
  2. Defines five callback functions that know how to update files and call the right  
  pipeline functions                                                                   
  3. Passes those callbacks to PseudoLabelingSetupWindow, which wires them to its      
  buttons     
  
   No major bugs. Two minor issues:                            
  1. Line 6410-6412: run_pseudolabeling doesn't return out_tsv — it returns None       
  implicitly (no return statement at the end). The fallback out_tsv =                  
  files.get("pseudolabel_csv") at line 6412 always fires. Not a bug since it works, but
   the comment "if wrapper stores it instead of returning" suggests this was
  unintentional. You could add return outpath at the end of run_pseudolabeling for
  clarity, or just accept the current pattern.
  2. Line 6458: sample_method_folder — same scoping issue we fixed. After your nonlocal
   fix, this will correctly read the user-set folder. Before the fix, it was always    
  None here too, so CGA methods always fell back to the CSV directory.
                                                                                       
  Otherwise the function is well-structured — clean separation of callbacks, proper    
  runlog entries, Score B enrichment integrated as a post-processing step, and method
  v1 auto-save at the end.                                                             
                                                            
❯ About issue 1, can you point out where to return? I don't think run_pseudolabeling   
  capables of returning file path... I might me wrong. If you mean making              
  run_pseudolabeling supports returning outpath, I think current workout (it generates 
   tsv, we capture out_tsv then)                                                       
                                                                                       
⏺ You're right — the current workaround is fine. run_pseudolabeling stores             
  files["pseudolabel_csv"] = outpath at line 5230, and on_start reads it back at line
  6412. It's a two-step pass through the files dict instead of a direct return, but    
  both paths end up with the same value.
                                                                                       
  Adding a return outpath at the end of run_pseudolabeling (after line 5255) would be  
  cleaner, but it's not necessary — the current pattern works and changing it risks
  breaking the subwin.after(0, refresh_tree) scheduling at line 5255. Leave it as-is. 
  

# open_pl_to_trainable_modal
## line 5873 - 6042
@Henry: It's CGA to trainable window. Fixed one edge case bug by Claude Code

# try_pl_to_trainable
## line 6043 - 6050
@Claude Code:
try_pl_to_trainable resolves the selection via _get_selected_context(), grabs the    
  files dict, and passes everything to open_pl_to_trainable_modal.
  
# try_pl_to_trainable_old
## removed
@Henry: old manner climbing node to get sample information


# Prepare dataset window button/UI design
## line 6054 - 6108 
@Henry: fixed ttk inconsistent button style, and note one issue found by Claude Code.
@Claude Code:
  2. Line 6722: experiment_status_labels["Unassigned"] = (exp_status_var,
  sample_status_var) — hardcoded to "Unassigned". When a user loads a named experiment,
   the status labels are never registered under that experiment name. That's why
  update_status_display does nothing for real experiments — the lookup at
  experiment_status_labels.get(exp_name) returns None.

  2. The fix would be to register under the current experiment name dynamically, or
  simply register under a generic key and have update_status_display always use that
  key. But this is tied to the broader update_status_display limitation we noted
  earlier — mark it for the future status bar redesign.
  3. Mixed tk.Button and ttk.Button: Lines 6706 and 6732-6733 use ttk.Button while
  everything else uses tk.Button. Not a bug but visually inconsistent — ttk buttons
  have different styling (themed/native look) vs tk buttons (classic look). Pick one
  for the whole panel.
 
# ML analysis Window 


# About ML method compatibility
@Claude Code:
  The clean approach when you add models: make BUILTIN_ML contain per-model-type       
  defaults:
                                                                                       
  BUILTIN_ML = {                                            
      "model": {
          "type": "RandomForest",                                                      
          "RandomForest": {"n_estimators": 400, "max_depth": None, ...},
          "XGBoost": {"n_estimators": 200, "learning_rate": 0.1, ...},                 
      },                                                                               
      "train": {...}
  }     
  
# _ml_state (dict initialization)
## line 6125 - 6150
@Henry: after first cleanup I noticed this is the ["editor"] printed left in param editor.

# _render_effective_from_state
## removed
@Claude Code: It was intended to update the ML Params editor's right-side
   preview pane (_ml_state["preview_txt"]) with the "effective" merged parameters.
   _emit_ml_state_changed_ui_refresh (which IS actively called) does the same job but
  better — it writes both editor and preview panes from _ml_state["editor"] and
  _ml_state["effective"].
 
# BUILTIN_ML (dict initialization)
## line 6155 - 6163
@Henry: removed redundant ones (there are 3 more) after fine-monitoring, moving to top of ml window

# _ensure_ml_state
## line 6165 - 6188
@Henry: As a fallback of if _ml_state is missing (first call or somehow getting dropped)
Debug lines added to monitor what's firing
 
 
# _write_editor_from_state
## removed
@Claude Code: It does the same job as  
  _emit_ml_state_changed_ui_refresh (line 6828): write _ml_state["editor"] and
  _ml_state["effective"] into the editor/preview Text widgets.                         
                                                                                       
  The difference: _write_editor_from_state uses _get_editor_preview_widgets() to find  
  the Text widgets, while _emit_ml_state_changed_ui_refresh reads them directly from   
  _ml_state["editor_txt"] and _ml_state["preview_txt"]. Same result, different lookup  
  mechanism.                                                

  _emit_ml_state_changed_ui_refresh is the active replacement. Delete                  
  _write_editor_from_state and also _get_editor_preview_widgets which is only called
  from this dead function.  
  
# _emit_ml_state_changed_ui_refresh
## line 6191 - 6211
@Claude Code: It writes both editor and preview panes from _ml_state["editor"] and
  _ml_state["effective"]
@Henry: asked for refactoring (removing redundant function covering this work), and changed the code
 t2.insert("1.0", json.dumps(snapshot_train_vars(), indent=2))  # 20260418, now right window read directly
so the right panel reads value from tk directly, not via _ml_state["editor"] (only lives in param editor)
  
# _write_train_ui_into_state
## removed
@Claude Code: (after refactoring) now _write_train_ui_into_state call becomes      
  pointless since nothing reads _ml_state["editor"] for training anymore. And the
  "Apply to ML Param[if not reflect]" button (line 8303) just triggers a preview       
  refresh, which is all it needs to do.

# _effective_params
## removed
@Henry: probably called by a prototype of _emit_ml_state_changed_ui_refresh

# _refresh_json_textboxes
## removed
@Henry: calls _effective_params and is confirmed a prototype of _emit_ml_state_changed_ui_refresh by Claude Code

# _deep_merge
## removed twice
@Claude Code:
  Both _deep_merge definitions are dead. You already have _merge_ml (line ~7823) which 
  does the same job with a cleaner variadic interface and is actively used throughout.
  Delete both _deep_merge definitions and merge_ml_params along with them.
  
# merge_ml_params
## removed
@Henry: deleted with second _deep_merge

# _merge_ml
## 3 candidates, only keep 1 and move to line 
    def _merge_ml(builtins: dict, exp_defaults: dict|None, editor: dict|None) -> dict:
        out = copy.deepcopy(builtins)
        for layer in (exp_defaults or {}, editor or {}):
            for k,v in layer.items():
                if isinstance(v, dict) and isinstance(out.get(k), dict):
                    out[k].update(v)
                else:
                    out[k] = v
        return out
    def _merge_ml(*layers):
        out = {}
        for layer in layers:
            if not layer: continue
            for k, v in layer.items():
                if isinstance(v, dict) and isinstance(out.get(k), dict):
                    out[k] = {**out[k], **v}
                else:
                    out[k] = v
        return out 
@Henry: analyzed result from Claude Code:
⏺ Three definitions:                                                                   
             
  ┌──────┬──────────────────────┬──────────────────────┬───────────────────────────┐   
  │ Line │      Signature       │     Merge style      │        Deep copy?         │   
  ├──────┼──────────────────────┼──────────────────────┼───────────────────────────┤   
  │      │ (builtins,           │                      │                           │   
  │ 7465 │ exp_defaults,        │ out[k].update(v)     │ Yes                       │   
  │      │ editor) — typed      │ (mutates in place)   │ (copy.deepcopy(builtins)) │   
  │      │ 3-arg                │                      │                           │
  ├──────┼──────────────────────┼──────────────────────┼───────────────────────────┤   
  │ 7735 │ (*layers) — variadic │ {**out[k], **v}      │ No (starts from {})       │
  │      │                      │ (creates new dict)   │                           │   
  ├──────┼──────────────────────┼──────────────────────┼───────────────────────────┤
  │ 7823 │ (*layers) — variadic │ {**out[k], **v}      │ No (starts from {})       │   
  │      │                      │ (creates new dict)   │                           │
  └──────┴──────────────────────┴──────────────────────┴───────────────────────────┘
  Lines 7735 and 7823 are identical. The second shadows the first. Both do shallow     
  one-level dict merge — nested dicts get their top-level keys merged but deeper
  nesting is replaced, not merged.                                                     
                                                            
  Line 7465 is the typed 3-arg version — functionally equivalent but limited to exactly
   3 layers. It deep-copies the first argument to avoid mutating the original, which is
   safer but unnecessary now since BUILTIN_ML is only read at init time. 
@Claude Code: keep the last one and move to top of ML window,


# Utility block clean up (probably by another ChatGPT session re-creating wheels)
  Delete all of these — zero callers, fully superseded:                                
                                                                                       
  ┌───────────────────────────────┬──────┬──────────────────────────────────────────┐  
  │           Function            │ Line │                 Why dead                 │  
  ├───────────────────────────────┼──────┼──────────────────────────────────────────┤  
  │ BUILTIN_ML (second def)       │ 7285 │ One of the extra definitions you're      │
  │                               │      │ deleting                                 │  
  ├───────────────────────────────┼──────┼──────────────────────────────────────────┤  
  │ _safe_read_json               │ 7294 │ Only called by the dead API functions    │  
  │                               │      │ below it                                 │  
  ├───────────────────────────────┼──────┼──────────────────────────────────────────┤  
  │ _safe_write_json              │ 7301 │ Same — only called by dead functions and │
  │                               │      │  snapshot_training_run                   │  
  ├───────────────────────────────┼──────┼──────────────────────────────────────────┤
  │ read_ml_defaults_from_exp     │ 7307 │ Zero callers — superseded by             │  
  │                               │      │ _read_ml_defaults_from_exp (line ~7446)  │
  ├───────────────────────────────┼──────┼──────────────────────────────────────────┤  
  │ write_ml_defaults_to_exp      │ 7312 │ Zero callers — superseded by             │  
  │                               │      │ _write_ml_defaults_to_exp                │
  ├───────────────────────────────┼──────┼──────────────────────────────────────────┤  
  │ read_ml_overrides_from_method │ 7321 │ Zero callers                             │
  ├───────────────────────────────┼──────┼──────────────────────────────────────────┤
  │ write_ml_overrides_to_method  │ 7326 │ Zero callers                             │
  ├───────────────────────────────┼──────┼──────────────────────────────────────────┤  
  │ load_ml_params_from_json      │ 7335 │ Zero callers — thin wrapper around       │
  │                               │      │ _safe_read_json                          │  
  ├───────────────────────────────┼──────┼──────────────────────────────────────────┤
  │ save_ml_params_to_json        │ 7338 │ Zero callers — thin wrapper around       │  
  │                               │      │ _safe_write_json                         │
  ├───────────────────────────────┼──────┼──────────────────────────────────────────┤
  │ _deep_merge (second def)      │ 7345 │ Already confirmed dead — _merge_ml       │
  │                               │      │ replaced it                              │  
  ├───────────────────────────────┼──────┼──────────────────────────────────────────┤
  │ merge_ml_params               │ 7355 │ Zero callers — superseded by _merge_ml   │  
  └───────────────────────────────┴──────┴──────────────────────────────────────────┘  
   
  Check before deleting:                                                               
                                                            
  ┌───────────────────────┬──────┬──────────────────────────────────────────────────┐  
  │       Function        │ Line │                     Verdict                      │
  ├───────────────────────┼──────┼──────────────────────────────────────────────────┤
  │ collect_versions      │ 7383 │ Check if snapshot_training_run is kept — if not, │
  │                       │      │  check other callers                             │
  ├───────────────────────┼──────┼──────────────────────────────────────────────────┤  
  │ snapshot_training_run │ 7408 │ Zero callers from train_model. Was designed for  │
  │                       │      │ provenance tracking but never wired in.          │  
  └───────────────────────┴──────┴──────────────────────────────────────────────────┘
                                                                                       
  collect_versions and snapshot_training_run are a nice provenance feature (save exact 
  environment + params + hashes alongside the model). But since snapshot_training_run
  was never called from train_model, they're dead. If you plan to add training         
  provenance tracking in the future, note them in your roadmap and delete — they're
  straightforward to recreate. If you want to keep them for near-term use, move them to
   the ML infrastructure block and wire snapshot_training_run into train_model after
  model save.

# _safe_read_json
## removed

# _safe_write_json
## removed
@Henry: absorbed into snapshot_training_run json saving part

# read_ml_defaults_from_exp
## removed

# write_ml_defaults_to_exp
## removed

# read_ml_overrides_from_method
## removed

# write_ml_overrides_to_method
## removed

# load_ml_params_from_json
## removed

# save_ml_params_to_json
## removed


# collect_versions
## line 6226 - 6250
@Henry: For collecting training environments (not parameters)

# snapshot_training_run
## line 6252 - 6277
@Henry: Saving report and enviroments to json file along with model file
Activated in 2026-04-18 code review; 
that's expected feature to allow diagnosis when sklearn throws version-dependent error
(happened once when I have all 3 devices - home PC/ lab PC/ lab Mac all in different pandas/sklearn version)
(perfect to mimic what will happen after we publish it -- WHY THE MODEL YOU PROVIDE SHOWS PREDICTION ERROR)


# Label consistency block
# _MANUAL_PURE
## line 6316
# _MANUAL_EXT
## line 6320

# _looks_tuplelike_for_structure
## line 6322 - 6328
@Henry: detect if it is tuple-like (need conversion, set to True)

# normalize_structure_for_training
## line 6330 - 6350
@Henry: Notice that it calls external module pretrain_normalizer.py. Need a check about string conversion then.

# load_model_any
## line 6355 - 6378
@Claude Code: A model loader that handles three file formats: .joblib, .pkl (both via         
  joblib.load), and .skops (via skops.io.load). Returns (model, loader_name).
@Henry: asked about skops fix, and now it should save skops as well to help avoid cross-version train error

# get_training_features
## line 6380 - 6400
@Claude Code: ecovers the exact list of feature column names that the model was trained on, avoid misalignment error

# balance_and_split
## line 6402 - 6554 
@Henry: Class Balance, Split, and Capping
@Claude Code:
  Step 0 — Drop rare classes (lines 7057-7061):                                        
  Remove any label with fewer than min_count samples from the entire dataset. This     
  prevents the splitter from failing when a class has too few samples to appear in both
   train and test. The dropped rows are recorded in info for transparency.             
                                                                                       
  Step 1 — Two modes diverge:                                                          
   
  Mode A: Cap before split (cap_training_only=False, lines 7088-7138):                 
  1. Find the largest minority class count (max_minor)      
  2. Compute cap: majority_factor × max_minor (e.g., if biggest minority has 100       
  samples and factor=3, cap Non-glycan at 300)                                  
  3. Randomly downsample the majority class to the cap                                 
  4. Split the already-capped data into train/test (and optionally val)
  5. Result: all three folds (train/val/test) reflect the capped distribution          
                                                                                       
  Mode B: Cap training only (cap_training_only=True, lines 7140-7193):                 
  1. Split uncapped data into train/test/val first                                     
  2. Then cap the majority class in the training fold only                             
  3. Test and val folds keep the natural (imbalanced) distribution                     
  4. Result: model trains on balanced data, but is evaluated on real-world-like        
  distribution                                                                         
                                                                                       
  The split mechanics (both modes):                                                    
  - First split: train_test_split(X, y, test_size=test+val) → produces train vs holdout
  - Second split: train_test_split(holdout, test_size=test/(test+val)) → splits holdout
   into val and test                                                                   
  - If val_size=0: skip second split, all holdout becomes test                         
  - stratify=y ensures class proportions are preserved in each fold; falls back to
  unstratified if stratification fails (too few samples per class)                     
                                                                                       
  Returns: (X_train, y_train, X_val, y_val, X_test, y_test, info) — six                
  DataFrames/Series plus a dict recording what happened (which mode, what was capped,  
  what was dropped). 
  
  
  
# _pull_state_into_vars
## removed
@Henry: for imprinting button state but  _render_effective_from_state gets removed so it is useless


# _guess_method_basename_for_pack
## line 6557 - 6620
@Claude Code:
⏺ We covered this during the initial audit. It's called from run_prediction to
  determine the output folder name for prediction results. The "guess" is because the
  prediction input CSV might not have an obvious link back to its method JSON.

  The function tries to find a method JSON associated with the prediction input so it
  can name the output folder meaningfully (e.g., sample_name.CGA.method/ instead of
  just predicted_output/). It searches in priority order:

  1. Column in the CSV containing a method path
  2. Exactly one *.method.json in the same folder as the input
  3. Multiple candidates → match by sample/experiment name from CSV columns
  4. Ask user via file dialog
  5. Last resort: use the CSV filename stem
  
# _pick_file_cli_or_gui
## line 6622 - 6638 
@Henry: Asked Claude Code. possible for CLI implementation in future.
Built 20250915. Called by create_unlabeled_from_method

# build_ml_params_panel
## removed
@Henry: migrated button to open_ml_params_window

# _get_current_context
## removed
@Henry: called by build_ml_params_panel, now useless 
(also it read exp json and set Title to None, probably a legacy-like behavior, since future ML params would be adapted to method json)

# _on_effective_params_ready
## removed
@Henry: called by build_ml_params_panel, now useless 


# _read_ml_defaults_from_exp
## line 6649 - 6655
@Henry: Read ML exp json
# _write_ml_defaults_to_exp
## line 6656 - 6667
@Henry: Write ML exp json
Note: what json is not clearly determined. Possible risk of erode standard 3-type json in GlycoMSP workflow
(when binding ML to dataset management, this part might need a clearer definition

# _debug_print_train_params
## line
@Henry: asked Claude Code to add debug prints since tkinter rendering bug is not resolved (impossible)
Note that Claude Code made mistake here 3 times (broken parethesis, missing value, one-liner)
Now I have adapted the print to 3 lines with end=""
Issue happened when context window is aroud 420k, so I believe 400k is a real potential bottleneck of window anxiety

# _get_editor_preview_widgets
## removed
@Henry: iirc the discussion with Claude Code said it monitor build_ml_params_panel's json fields


# _widget_alive -> widget_alive
## line 6720 - 6725
@ Henry: kept it when in future updates requires this functionality again. Migrating it to top so it can be called in all modules
@ Claude Code:
⏺ Correct. The two versions are identical logic — one is a class method (self, w), one 
  is a standalone function (w). Moving the standalone version to module-level as     
  widget_alive(w) makes it available to both MetadataEditorWindow and                  
  open_ml_analysis_window without duplication.
                                                                                       
  One detail: the MetadataEditorWindow version is a method called as                   
  self._widget_alive(w). After moving to module-level, those calls change to
  widget_alive(w) — drop the self. prefix. That's a minor refactor across ~5 call sites
   in the class.
   
# _update_ml_summary()  : 3 same functions, remove first 2 and refactor 3rd one
## line 6815 - 6830

@Claude Code: 
⏺ Three copies, the last one wins. Comparison:                                         
                                                                                       
  ┌────────────┬─────────────────────┬─────────────────────┬──────────────────────┐    
  │            │   Line 7674 (1st)   │   Line 7715 (2nd)   │   Line 7835 (3rd,    │    
  │            │                     │                     │       active)        │    
  ├────────────┼─────────────────────┼─────────────────────┼──────────────────────┤    
  │ Exception  │ ✅ try/except with  │ ❌ None             │ ❌ None              │    
  │ handling   │ fallback            │                     │                      │    
  ├────────────┼─────────────────────┼─────────────────────┼──────────────────────┤
  │ Source     │ effective_ml_params │ effective_ml_params │ effective_ml_params  │    
  │            │  only               │  or {}              │ or BUILTIN_ML        │    
  ├────────────┼─────────────────────┼─────────────────────┼──────────────────────┤
  │            │                     │                     │ model + train        │    
  │ Shows      │ model params only   │ model params only   │ (n_estimators,       │    
  │            │ (5 fields)          │ (5 fields)          │ class_weight,        │
  │            │                     │                     │ test_size, val_size) │    
  ├────────────┼─────────────────────┼─────────────────────┼──────────────────────┤
  │ BUILTIN_ML │ No — falls back to  │ No — empty dict if  │ Yes — reads          │
  │  fallback  │ hardcoded 400 in    │ no effective        │ BUILTIN_ML if no     │    
  │            │ .get()              │                     │ effective            │
  └────────────┴─────────────────────┴─────────────────────┴──────────────────────┘ 

# _on_change(*_)
## line 6192 - 6196
@Henry: moved to near by the top after confirming with Claude Code

# for v in () trace.add
## line 6197 - 6202
@Henry: noticed this is live update trace in open_ml_analysis_window, moved neap top after variables are initialized


# tk variable block
## line 6172 - 6188 
@Henry: defines all values for ML, shared on top (moved from middle and removed duplicates)


# open_ml_params_window
## line 6821 - 7082
@Henry: asked Claude Code to confirm duplicates inside and remove them
Absorbed build_ml_params_panel and adjusted the button and window size for better UX
(Large reconstruction and refactoring conducted, check session log would be better if looking up what's happened)

# open_train_settings
## line 7084 - 7229
@Henry: asked Claude Code to fix several issues, and now debug prints are fine.
Despite tk/ttk button rendering issue, now params are correctly stored, updated and printed as debug lines

# train model
## line 7231 - 7614
@Henry: one-line summary: fixed outputs and mass feature tracing. Support both full RF with w/o Non-glycan report
@Henry: asked stage-separated summary for code review 
@Claude Code:
⏺ Here's the stage summary:
                                   
  Stage: 1) Import & config                                                            
  Lines: 7337-7372                                                                     
  Operation: Import sklearn/joblib/normalizer. Read use_mass_feature. Resolve cfg from 
    effective_ml_params or snapshot_train_vars(). Debug prints.                        
  Key action: cfg = final ML parameters                                                
  ────────────────────────────────────────                                             
  Stage: 2) Load & validate                                                         
  Lines: 7377-7392                                                                     
  Operation: Check train_csv_path exists. Validate label_col. Load CSV.            
  Key action: df = training DataFrame                                                  
  ────────────────────────────────────────                                         
  Stage: 3) Normalize labels                                                           
  Lines: 7394-7419                                                                 
  Operation: If label is "Structure", run normalize_structure_for_training via         
    pretrain_normalizer. Convert to string series y_raw.                           
  Key action: Labels standardized to F1H5N4S2 format                                   
  ────────────────────────────────────────                  
  Stage: 4) Balance & split                                                            
  Lines: 7421-7451
  Operation: Call balance_and_split with or without balancing. Mode A (cap before      
  split)                                                    
    or Mode B (cap training only).
  Key action: X_train/val/test, y_train/val/test, info
  ────────────────────────────────────────
  Stage: 5) MS1 feature handling                                                       
  Lines: 7454-7484
  Operation: _prep_ms1_block: include/exclude protonatedmass. Non-glycan rows → 0.0.   
  NaN                                                       
    fill. Apply to all three splits.
  Key action: Feature matrix cleaned
  ────────────────────────────────────────
  Stage: 6) Label encoding                                                             
  Lines: 7486-7490
  Operation: LabelEncoder.fit_transform on train, .transform on val/test.              
  Key action: y_train_enc, y_val_enc, y_test_enc            
  ────────────────────────────────────────
  Stage: 7) Build & train model                                                        
  Lines: 7492-7526
  Operation: Build RandomForestClassifier from cfg["model"]. _drop_meta_and_get_X to   
    strip non-features. Align column order. model.fit().    
  Key action: Model trained
  ────────────────────────────────────────
  Stage: 8) Evaluate                                                                   
  Lines: 7529-7616
  Operation: Raw predictions on test. Binary glycan-vs-non report. Optional confidence 
    threshold demotion. Full classification report. Glycan-only report. Confusion
    matrix.
  Key action: Performance metrics
  ────────────────────────────────────────
  Stage: 9) Notify user                                                                
  Lines: 7618-7641
  Operation: Build summary message with cap info, thresholding status, full report.    
    messagebox.showinfo.                                    
  Key action: User sees results
  ────────────────────────────────────────
  Stage: 10) Save artifacts                                                            
  Lines: 7643-7688
  Operation: Save model (.skops or .joblib fallback), label encoder (.joblib),         
    performance report (.txt), feature list (.json), training run snapshot (.json).
  Key action: 5 files written
  ────────────────────────────────────────
  Stage: 11) Persist to exp                                                            
  Lines: 7691-7717
  Operation: Write train_parameters dict back to linked .exp.json if available.        
  Key action: Experiment JSON updated                       

⏺ _prep_ms1_block is called three times at lines 7481-7483:                            
                                                                                       
  X_train = _prep_ms1_block(X_train, y_train, use_mass_feature)                        
  X_val   = _prep_ms1_block(X_val,   y_val,   use_mass_feature)                        
  X_test  = _prep_ms1_block(X_test,  y_test,  use_mass_feature)                        
                                                                                       
  X — the feature DataFrame for one split (train, val, or test). Columns are ion mass  
  features (like "344.1726", "464.249") plus possibly "protonatedmass" and metadata    
  columns.                                                                             
                                                            
  y — the corresponding label Series (like "F1H5N4S2", "Non-glycan", etc.). Still      
  strings at this point — LabelEncoder runs after this step.
                                                                                       
  use_mass — boolean from the "Use protonated mass as a model feature" checkbox.       
   
  What it does:                                                                        
                                                            
  - If use_mass=True and protonatedmass exists: coerce to numeric, set Non-glycan rows 
  to 0.0 (because Non-glycan scans don't have meaningful precursor masses — using the
  real mass would leak information), fill remaining NaNs with 0.0.                     
  - If use_mass=True but no protonatedmass column: print warning, continue without it.
  - If use_mass=False: drop protonatedmass entirely so the model never sees it.        
  - Finally: fillna(0.0) on all features as a safety net.                              
                                                                                       
  Why y is needed: Line 7467-7468 uses y to identify which rows are Non-glycan to zero 
  out their mass. Without y, it couldn't distinguish glycan from non-glycan rows. The  
  label is still in string form here, so it checks y_ser.eq("Non-glycan") directly. 

  maj_idx = list(le.classes_).index(majority_label_var.get())                          
          y_true_bin = (y_test_enc != maj_idx).astype(int)   # 1=glycan, 0=non         
          y_pred_bin = (y_pred      != maj_idx).astype(int), what do these lines do?   
                                                                                       
⏺ This creates a simplified binary evaluation — "glycan vs non-glycan" regardless of   
  how many glycan compositions exist:                                                
                                                                                       
  maj_idx — finds the integer index of the majority label (usually "Non-glycan") in the
   encoder's class list. If le.classes_ is ["F1H5N4S2", "H5N4S2", "Non-glycan"], then  
  maj_idx = 2.                                                                         
                                                            
  y_true_bin — converts true test labels to binary: anything that's NOT the majority   
  class gets 1 (glycan), the majority class gets 0 (non-glycan). So all glycan
  compositions collapse into one group.                                                
                                                            
  y_pred_bin — same conversion on predictions.                                         
   
  Then lines 7541-7543 compute binary precision/recall/F1/accuracy on this simplified  
  view:                                                     
                                                                                       
  p,r,f,_ = precision_recall_fscore_support(y_true_bin, y_pred_bin, average="binary")
                                                                                       
  This answers: "ignoring which specific glycan it is, how well does the model         
  distinguish glycan from non-glycan?" — a useful high-level metric alongside the      
  per-class report. A model might misclassify F1H5N4S2 as H5N4S2 (wrong composition but
   still glycan) — the per-class report penalizes this, but the binary report counts it
   as a correct glycan detection.


# A logical/design issue about params passing to train/split/val
@Henry: Got it. So this is benifited from completence of snapshot_train_vars() having all    
  information we need? Is the eff sync with tk variables? I mean, won't directly read  
  the tk variables bypass the if effective_ml_params there? Although I presumably      
  think the values are synced (not sure, really) so it should be fine...               
                                                                                       
@Claude Code: You're right to be concerned. They're not always synced. Here's the gap:                   
  - tk variables — updated by sliders/checkboxes in real-time                          
  - effective_ml_params — only updated when user clicks "Validate & Use" in ML Params
  editor                                                                               
                                                            
  If the user opens the ML Params editor, edits the JSON to set n_estimators=200,      
  clicks Validate & Use → effective_ml_params has n_estimators=200. Then they go to
  Train Settings, drag the n_estimators spinbox to 800 → tk variable is now 800, but   
  effective_ml_params still has 200.                        

  When train_model runs:
  - cfg = effective_ml_params (has n_estimators=200)
  - rf = cfg["model"] → n_estimators=200 for the model                                 
  - But balance_and_split reads test_split_var.get() directly from tk → gets the slider
   value                                                                               
                                                                                       
  So the model params come from the editor, but the split/balance params come from the 
  sliders. They're two different sources in the same training run.                     
                                                            
  This is actually fine for your current design — the editor is for advanced model     
  tuning, the sliders are for split/balance tuning. But it's important to know they can
   diverge. If you want full consistency, train_model should read everything from cfg  
  instead of tk variables. But that's a bigger refactor — note it for now.
  
  

⏺ The thresholding block is a post-prediction filter that demotes low-confidence glycan
   predictions to Non-glycan. Step by step:                                          
                                                                                       
  tau — confidence threshold (e.g., 0.65). If the model's best class probability is    
  below this, the prediction is considered uncertain.                                  
                                                                                       
  margin — how close the Non-glycan probability must be to the best probability (e.g., 
  0.05). This checks if the model is "nearly indifferent" between its top pick and
  Non-glycan.                                                                          
                                                            
  The demotion logic (line 7568):
  demote = (y_pred != majority_idx) & (maxp < tau) & (p_major >= (maxp - margin))
                                                                                 
  Three conditions must ALL be true to demote:                                         
  1. The prediction is a glycan (not already Non-glycan)                               
  2. The model's confidence is below tau (uncertain)                                   
  3. The Non-glycan probability is within margin of the best probability (model almost 
  picked Non-glycan anyway)                                                            
                                                                                       
  If all three pass, the prediction gets overwritten to Non-glycan.
                                                                                       
  Why it worsens results: It's designed to reduce false positive glycan predictions —  
  but with small datasets, the model's probability estimates are poorly calibrated.    
  Many correct glycan predictions have low confidence (the model is right but unsure), 
  so the threshold demotes true positives alongside false positives. The net effect:
  fewer false positives but also fewer true positives, lowering recall and often F1.

  It's more useful with larger, well-calibrated models. With your current dataset size,
   keeping it off is the right call.
   
   
❯ Thanks for the call. GPT strongly recommended me to add this even I explained we're  
  on a small datasets then. I'm happy that disabling it is on the right track. Now can 
   you tell me what these lines do:        used = unique_labels(y_test_enc, y_pred)    
          used_names = [le.classes_[i] for i in used]                                  
                                                                                       
⏺ unique_labels(y_test_enc, y_pred) — finds all class indices that appear in either the
   true labels OR the predictions. Some classes might exist in training but never show 
  up in the test set or predictions (too rare, or thresholding removed them). This     
  collects only the classes that actually participate in the evaluation.
                                                                                       
  used_names — converts those integer indices back to human-readable labels using the  
  encoder. So [0, 2, 5] becomes ["F1H5N4S2", "H5N4S2", "Non-glycan"].
                                                                                       
  These are passed to classification_report(y_test_enc, y_pred,                        
  target_names=used_names) so the report only shows rows for classes that actually
  appear — avoids cluttering the report with classes that have 0 samples in both true  
  and predicted, which would show 0.00 across all metrics.  

❯ Got it. Then what do these following lines do?         labels_all =                  
  list(le.classes_)                                                                    
          cm = confusion_matrix(y_test_enc, y_pred, labels=np.arange(len(labels_all))) 
          no_pred = [labels_all[j] for j, s in enumerate(cm.sum(axis=0)) if s == 0]    
          no_true = [labels_all[i] for i, s in enumerate(cm.sum(axis=1)) if s == 0]    
                                                                                       
⏺ This checks for "invisible" classes — ones that exist in training but vanished from
  test or predictions:
                                                                                       
  cm = confusion_matrix(...) — builds the full confusion matrix including ALL classes  
  (not just the ones in used). Each row = true class, each column = predicted class.   
                                                                                       
  no_pred — classes where the column sums to 0: the model never predicted this class   
  for any test sample. The class exists in training but the model never outputs it —
  could mean the class is too similar to another and always gets confused.             
                                                            
  no_true — classes where the row sums to 0: no test samples have this true label. The 
  class was in training but the random split put zero samples in the test set — happens
   with rare classes.                                                                  
                                                            
  Both are printed as warnings (lines 7593-7594) so you know which classes are missing 
  from the evaluation. If no_true lists a class, the classification report can't
  evaluate it (no ground truth). If no_pred lists a class, the model has a blind spot  
  for it.                                                   

❯ Oh, that generates the debug lines I'm seeing, right? [ML] No predicted samples for: 
   ['F1H5N3', 'F1H5N3S1G2', 'F1H5N6G2', 'F1H6N2', 'F1H6N6', 'F2H5N4', 'F3H4N5G1',      
  'H4N5S1', 'H5N3', 'H8N6', 'H9N5']                                                    
  [ML] No true samples in test for: []                                                 
                                                                                       
⏺ Yes, exactly. The model knows those 11 compositions from training but never predicts 
  any of them on test data — they always get classified as something else (likely a  
  similar composition or Non-glycan). No true samples missing from test means the split
   was good — every class has at least some test representation.
                                                                                       
  Those 11 "blind spot" classes are likely rare compositions with very few training    
  samples. The model doesn't have enough examples to learn their patterns confidently,
  so it defaults to a more common class. This is expected with small datasets and many 
  classes — not a bug, just a data limitation. More training samples for those
  compositions would fix it.

⏺ This block produces a glycan-only classification report — evaluating only how well
  the model distinguishes between different glycan compositions, excluding Non-glycan
  entirely:                                                                            
   
  Lines 7599-7601: Find the majority class index (Non-glycan). try/except in case the  
  majority label somehow isn't in the classes.              
                                                                                       
  Line 7604: mask_glycan_true = (y_test_enc != maj_idx) — boolean mask selecting only  
  test samples whose TRUE label is a glycan (not Non-glycan).
                                                                                       
  Lines 7606-7607: Build the list of glycan-only class indices and names — everything  
  except the majority class.
                                                                                       
  Lines 7608-7614: classification_report on only the glycan-masked subset:             
  classification_report(
      y_test_enc[mask_glycan_true],   # true labels, glycan rows only                  
      y_pred[mask_glycan_true],       # predicted labels for those same rows
      labels=glycan_label_indices,     # only evaluate glycan classes                  
      target_names=glycan_names,       # human-readable names                          
      zero_division=0                                                                  
  )                                                                                    
                                                                                       
  Why this matters: The full report includes Non-glycan, which is usually the largest  
  class with high precision/recall — it inflates the overall accuracy. The glycan-only 
  report shows how well the model discriminates between specific compositions (e.g.,
  does it confuse F1H5N4S2 with H5N4S2?). This is the scientifically meaningful metric 
  — glycan identification accuracy, not just glycan detection.

  Line 7615: Printed to terminal only, not included in the messagebox or saved report  
  file. You might want to save it alongside the main report if it's useful for the
  manuscript.    
  
  balance_and_split returns different key names depending on which mode ran:           
  - Mode A (cap before split): info["majority_cap_applied_to"] and info["majority_cap"]
  - Mode B (cap training only): info["train_majority_cap_applied_to"] and              
  info["train_majority_cap"]                                             
                                                                                       
  The or chains check both key names so the message works regardless of mode. Then
  cap_prefix distinguishes "Capped" (all data) vs "Capped training" (train fold only). 
  The result is a single human-readable line like:
                                                                                       
  - "Capped training 'Non-glycan' at 450" (Mode B)                                     
  - "Capped 'Non-glycan' at 300" (Mode A)
  - "No majority cap applied." (balancing disabled)                                    
                                                            
  This goes into the msg list for the completion messagebox. 

# train_with_optional_ng
## removed
@Henry: previous code for injecting negatives when easy-negatives/Non-glycans hasn't been introduced

# functions for prediction

# select_model_file
## line 7617 - 7622
@Henry: function for seleting model file
# select_predict_input
## line 7624 - 7629
@Henry: function for selecting predicable csv
@Henry: consider both moving up to near top of open_ml_analysis_window, that's generic utility



# _read_any_table
## removed 
@Henry: another duplicate before run_prediction

# _canon_comp
## removed
@Henry: the real active one is in run_prediction, delete this one in open_ml_analysis_window

# _find_col, _LABEL_RE, _canon_label, _label_from_counts, apply_mass_gate.  
## removed  
@Claude Code: Keep _canon_comp and _canon_safe — they're the live versions.

# _canon_comp +  _comp_pat definition
## line 7789 - 7805 
@Henry: accumulates str->int compositional count to dict, and return FxHxNxSxGxKDNxAxsxpx string format

# _clean_cols
## moved to top utility and removed one in prepare dataset window
@Claude Code:It strips whitespace and removes the Unicode BOM character (\ufeff) from all column headers
@Henry: yeah I remembered this was a fix to something invisible (from human)
@Claude Code: promote it to place near robust_read_csv


# run_prediction
## line  7637 - 8245
@Henry: Too long, see doc highlights below

⏺ Here's the stage summary:                                                            
                                                                                       
  ┌───────────────┬───────────┬────────────────────────────────────────────────────┐   
  │    Stage     │   Lines   │                      Operation                      │   
  ├──────────────┼───────────┼─────────────────────────────────────────────────────┤   
  │ 1) Import &  │           │ Import joblib/numpy. Check model and input paths    │   
  │ validate     │ 7790-7826 │ exist. Load model via load_model_any. Load CSV.     │   
  │              │           │ Load label encoder. Handle version mismatch errors. │   
  ├──────────────┼───────────┼─────────────────────────────────────────────────────┤   
  │              │           │ Get training feature list. Add placeholder columns  │   
  │ 2) Feature   │ 7829-7861 │ if model expects MS1/delta_ppm but input lacks      │   
  │ alignment    │           │ them. Build X, align to training feature order,     │
  │              │           │ fill missing with 0.                                │   
  ├──────────────┼───────────┼─────────────────────────────────────────────────────┤
  │ 3) Predict   │ 7862-7863 │ model.predict(X) — raw predictions.                 │
  ├──────────────┼───────────┼─────────────────────────────────────────────────────┤   
  │ 4)           │           │ Optional margin-aware threshold: demote             │
  │ Confidence   │ 7870-7908 │ low-confidence glycan predictions to Non-glycan.    │   
  │ demotion     │           │ Note: always reads tau/margin regardless of         │
  │              │           │ enable_thresh_var — no enable check here.           │   
  ├──────────────┼───────────┼─────────────────────────────────────────────────────┤
  │ 5) Label     │           │ Convert integer predictions back to string labels   │
  │ decode       │ 7911-7920 │ via LabelEncoder or model.classes_. Assign to       │   
  │              │           │ df['Predicted_Label'].                              │
  ├──────────────┼───────────┼─────────────────────────────────────────────────────┤   
  │ 6) Gate      │           │ Ensure ppm_precursor and pred_ok columns exist with │
  │ column       │ 7922-7928 │  defaults (NaN and True).                           │   
  │ defaults     │           │                                                     │
  ├──────────────┼───────────┼─────────────────────────────────────────────────────┤   
  │ 7) Reporter  │ 7930-7937 │ Collect class names, compute probabilities          │
  │ prep         │           │ (duplicate of earlier block).                       │   
  ├──────────────┼───────────┼─────────────────────────────────────────────────────┤
  │ 8)           │           │ Define _canon_comp, _find_col, _canon_label,        │   
  │ Composition  │ 7942-8099 │ _label_from_counts, apply_mass_gate — all for       │   
  │ helpers      │           │ precursor gate matching. apply_mass_gate defined    │
  │              │           │ but never called.                                   │   
  ├──────────────┼───────────┼─────────────────────────────────────────────────────┤
  │ 9) Packaging │           │ Determine MS1/gate suffixes.                        │
  │  setup       │ 8134-8149 │ _guess_method_basename_for_pack for output folder.  │
  │              │           │ Create output directory.                            │   
  ├──────────────┼───────────┼─────────────────────────────────────────────────────┤
  │ 10)          │           │ If gate enabled: load in-silico CSV, build          │   
  │ Precursor    │ 8151-8539 │ composition→mass lookup, compute PPM for each       │   
  │ gate         │           │ glycan prediction, mark pass/fail, write debug      │
  │              │           │ CSVs, gate counts.                                  │   
  ├──────────────┼───────────┼─────────────────────────────────────────────────────┤
  │ 11)          │           │ Normalize proba_vector lengths across rows for the  │
  │ Probability  │ 8542-8570 │ reporter.                                           │   
  │ harmonize    │           │                                                     │
  ├──────────────┼───────────┼─────────────────────────────────────────────────────┤   
  │ 12) Report   │ 8572-8604 │ summarize_predictions + write_prediction_report +   │
  │              │           │ show_prediction_summary_popup.                      │   
  ├──────────────┼───────────┼─────────────────────────────────────────────────────┤
  │ 13) Save &   │           │ Save predicted CSV. Three duplicate                 │   
  │ notify       │ 8593-8617 │ messagebox.showinfo calls — lines 8607, 8608, and   │   
  │              │           │ 8617.                                               │
  └──────────────┴───────────┴─────────────────────────────────────────────────────┘   
   


⏺ Step by step:                                                                        
                                         
  X = df.drop(columns=["MS2scan_no"], errors="ignore") — Start building the feature    
  matrix from the prediction CSV. Drop the scan ID column (not a feature). Everything
  else stays for now.                                                                  
                                                            
  missing = sorted(set(train_feats) - set(X.columns)) — Find columns the model expects 
  but the input CSV doesn't have. Example: model trained with 280 ion masses, but this
  CSV only has 275 — the 5 missing ones go in this list.                               
                                                            
  extra = sorted(set(X.columns) - set(train_feats)) — Find columns the input CSV has   
  but the model doesn't expect. Example: the CSV has "Structure", "composition" — text
  columns the model never trained on.                                                  
                                                            
  for c in missing: X[c] = 0.0 — Add the missing columns filled with zeros. The model  
  needs exactly the same columns as training — a missing ion mass feature means "this
  ion was never observed" which is zero intensity.                                     
                                                            
  X = X[list(train_feats)] — Reorder columns to match the exact training order AND drop
   extra columns in one step. Only keeps columns in train_feats, in that order. This is
   critical because RandomForest expects features in the same position as during       
  training — column names alone aren't enough, the order must match.

  So the whole block ensures: same columns, same order, no extras, no missing — the    
  prediction input is exactly shaped like the training input.
                                                                                       
❯ Thanks. So to prevent errors we reordered the column order if needed (this is fine), 
   and the missing feature values in prediction data gets 0 -- I think I should log    
  this, since that reduces the ability of model to predict properly (Not model's       
  fault, it's mistake from input data). A question out of the scale here: how do       
  prediction results get the proper tracable MS2scan_no and probably unique id or      
  something aligning back?                                                             
                                                                                       
⏺ The df (original DataFrame) is never modified by the feature alignment — only X (the 
  copy) gets columns dropped/reordered. After y_pred = model.predict(X), the
  predictions are assigned back to the original: df['Predicted_Label'] = y_pred (line  
  7920). Since X was built from df without changing row order or index, row 0 in y_pred
   corresponds to row 0 in df, which still has MS2scan_no and all other original       
  columns.                                                  

  So the alignment is by row position — df keeps all identity columns intact, X is a   
  stripped-down copy used only for prediction, predictions map back by index. The final
   saved CSV (df.to_csv(...)) has both MS2scan_no and Predicted_Label side by side.    
                                                            
❯ Ok, so can I recognize it as: we slice the df with fixed order (no index reset, no   
  sorting), and only feed model with information it knows, and after model tells us    
  the prediction result, we concate the result back to original df with sliced df      
  having prediction label. Something like that?                                        
                                                                                       
⏺ Exactly. The flow is:                                                                
                                                                                       
  1. df = full CSV (scan IDs, features, metadata — everything)                         
  2. X = slice of df (features only, reordered to match training) — same row order, no
  index reset                                                                          
  3. y_pred = model.predict(X) — one prediction per row, same order
  4. df['Predicted_Label'] = y_pred — stick predictions back onto the full DataFrame   
                                                            
  Row 0 in X came from row 0 in df, so prediction 0 belongs to row 0 in df. No join or 
  merge needed — just positional alignment. The original MS2scan_no, protonatedmass,
  and everything else in df are untouched throughout.
  
⏺ if le is not None: — LabelEncoder was loaded from the _labelencoder.joblib sidecar   
  file. Use le.inverse_transform to convert integer predictions back to string labels
  (e.g., 2 → "F1H5N4S2"). The bare except: pass silently ignores errors — if the       
  demotion block already converted predictions to strings via le.inverse_transform at  
  line 7906, this call would try to double-decode and fail, which the pass catches.    
                                                                                       
  elif hasattr(model, 'classes_'): — No label encoder file found, but the model itself 
  stores class names (sklearn models do this). Map integers to strings using
  model.classes_[i]. The isinstance(i, int) check handles the case where predictions   
  are already strings (passes them through).       

  Takes a DataFrame and a list of candidate column names. Builds a lowercase lookup of 
  all columns, then checks each candidate. Returns the original column name (preserving
   case) for the first match. Returns None if no candidate matches.                    
                                                            
  Example: _find_col(df, ["Predicted","Predicted_Label","label"]) — if the DataFrame   
  has "Predicted_Label", returns "Predicted_Label". Handles CSVs from different sources
   that might capitalize headers differently.                                          
                                                            
  Same pattern as the aliases dict and has_lower we saw in                             
  build_trainable_from_pseudolabels, just packaged as a reusable function.
  
⏺ _LABEL_RE is a regex that matches composition labels in strict FHNSGKDN order using
  named groups. Unlike _comp_pat which finds tokens in any order, this requires them in
   sequence:                                                                           
   
  - (?:F(?P<F>\d+))? — optional Fuc with count, captured as group F                    
  - (?:H(?P<H>\d+))? — optional Hex, captured as H          
  - ...and so on through KDN                                                           
                                                            
  Each ? makes the group optional, so "H5N4" matches (no F, S, G, KDN).                
                                                            
  _canon_label uses it:                                                                
                                                            
  m = _LABEL_RE.fullmatch(s)                                                           
  F = int(m.group("F") or 0)                                
  H = int(m.group("H") or 0)
  # ... extract all counts ...                                                         
  # rebuild in canonical order, omitting zeros
                                                                                       
  Difference from _canon_comp:                              
  - _canon_comp uses findall — tolerates any order, handles modifiers (A/s/p)          
  - _canon_label uses fullmatch — requires strict FHNSGKDN order, no modifiers, fails  
  on unrecognized formats                                                            
                                                                                       
  _canon_label is stricter because it's used for the precursor gate where prediction
  labels (always output in canonical order from the model) need to match library       
  entries. _canon_comp is more forgiving for parsing messy input from various sources.
                                                                                       
  In practice they produce the same output for well-formed labels — _canon_label is    
  just faster to reject garbage strings since fullmatch fails immediately rather than
  scanning for tokens.  
  
# 1) class names
⏺ Gets the list of all class names the model knows (e.g., ["F1H5N4S2", "H5N4",         
  "Non-glycan", ...]). Prefers the label encoder's le.classes_ (ordered list matching
  integer indices), falls back to model.classes_ if no encoder. This list is needed    
  later by the prediction reporter to label probability columns and generate per-class 
  summaries.   
  
proba_cols = None 
⏺ It checks if probability data exists and is a usable array-like object:              
                                         
  - proba is not None — model has predict_proba and it returned something              
  - isinstance(proba, (list, tuple)) — it's a Python list/tuple (unlikely but
  defensive)                                                                           
  - hasattr(proba, "shape") — it's a numpy array (the normal case from sklearn)
                                                                                       
  If true, it tries to attach probabilities to pred_df in one of two formats:          
                                                                                       
  Primary (line 8124): Store each row's probability vector as a single cell            
  (proba_vector column). This is compact — one column holds a list per row like [0.05, 
  0.85, 0.10].                                                                         
                                                            
  Fallback (lines 8126-8131): If that fails, expand into separate columns — one per    
  class like proba_F1H5N4S2, proba_Non-glycan, etc. More verbose but easier to filter
  in Excel.                                                                            
                                                            
  The reporter later uses whichever format was produced to compute confidence metrics, 
  entropy, margin, etc. In practice, with sklearn's RF, proba is always a numpy array
  so the hasattr(proba, "shape") path fires and the primary format succeeds.     
  
 ❯ okay. Now digging a bit on the column names, does the space in string also get       
  stripped? such as "theoretical mass".                                                
                                                                                       
⏺ No — .strip() only removes leading/trailing whitespace, not internal spaces.         
  "theoretical mass" stays as "theoretical mass", and " theoretical_mass " becomes     
  "theoretical_mass".                                                                  
   
  But that's not a problem here — the lookup candidates use underscores                
  ("theoretical_mass", "mass", "composition"), and CSV column headers from compnewv4
  use underscores too. If someone manually created a CSV with "theoretical mass" (space
   instead of underscore), the lookup would miss it.        

  That's an edge case not worth handling now — your pipeline always produces           
  underscore-separated headers. If it ever becomes an issue, add c.replace(" ", "_") to
   clean_columns.   
   
 comp2masses
 ⏺ Builds a dictionary mapping each composition label to its list of theoretical masses:
            
  - Iterates through the library rows, pairing comp_col (composition string) with      
  mass_col (theoretical mass)                                                          
  - Skips rows where mass is NaN                                                       
  - Canonicalizes the composition via _canon_safe (ensures matching format with        
  prediction labels)                                                                   
  - setdefault(comp_key, []).append(float(mass)) — if this composition hasn't been     
  seen, create an empty list; then append the mass                                     
                                                            
  Result example:                                                                      
  {                                                         
      "F1H5N4S2": [2076.0234, 2076.0235],  # same composition, slightly different 
  masses (different adducts/modifications)                                        
      "H5N4": [1663.8123],                                                             
      "F1H5N4S1": [1785.9012],
  }                                                                                    
                                                                                       
  A composition can have multiple masses because different modification combinations   
  (sulphate positions, linkage variants) produce the same FHNSGKDN label but slightly  
  different theoretical masses. The gate later checks if the observed precursor mass is
   within PPM tolerance of any of the listed masses for that composition.   
   
⏺ Correct. The logical order should be:

  1. Initialize gate columns (ppm_precursor = NaN, pred_ok = False) — already in place 
  2. Compute PPM for glycan rows (the __canon_pred__ → loop → best_ppm → keep block) —
  this is where real values get assigned                                               
  3. Write back results (pred_df.loc[idx, "ppm_precursor"] = best_ppm, pred_df.loc[idx,
   "pred_ok"] = keep)                                                                  
  4. Then run GATE DIAGNOSTICS (now kept/dropped reflect reality)
  5. Then run EXTRA DIAGNOSTICS (no_lut rows with counts-match)                        
  6. Save gated output CSV and counts                                                  

⏺ Preparing the arrays for the PPM computation loop:

  - pred_df["__canon_pred__"] — Canonicalizes every prediction label via _canon_safe   
  (e.g., "H5F1N4S2" → "F1H5N4S2"). Stored as a temp column for matching against
  comp2masses.                                                                         
  - gly_mask = (~is_ng).to_numpy() — Boolean array: True for glycan rows, False for
  Non-glycan. Inverts is_ng we created earlier. .to_numpy() converts from pandas Series
   to numpy array for faster indexing.
  - idx = np.nonzero(gly_mask)[0] — Integer indices of all glycan rows (e.g., [0, 1, 3,
   7, ...]). Used later to write results back to the correct positions in pred_df via  
  pred_df.loc[idx, ...].
  - labels — Canonical composition strings for glycan rows only (numpy array of        
  strings).                                                                            
  - obs — Observed precursor masses for glycan rows only (numpy array of floats).
  errors="coerce" converts non-numeric values to NaN.                                  
  - best_ppm — Pre-allocated array of NaN, one per glycan row. Will hold the best
  (closest) PPM error found for each row during the loop.                              
  - keep — Pre-allocated array of False, one per glycan row. Will be set to True for
  rows that pass the PPM threshold.                                                    
                                                            
  So: extract glycan-only data into parallel arrays, pre-allocate results arrays, then 
  the loop fills them in. After the loop, pred_df.loc[idx, "ppm_precursor"] = best_ppm
  and pred_df.loc[idx, "pred_ok"] = keep write results back to the correct rows.
  
⏺ _parse_counts and _masses_by_counts serve a different purpose than _canon_comp:

  _parse_counts(label) — Takes a label string like "F1H5N4S2s1" and returns a dict of
  individual counts: {"F":1, "H":5, "N":4, "S":2, "G":0, "K":0, "A":0, "s":1, "p":0}.
  It breaks the string into components — not rebuilding a canonical string, just
  extracting numbers.

  _masses_by_counts(lbl) — Takes a label, parses it to counts via _parse_counts, then
  searches lib_counts for rows where every count column matches exactly. Returns the
  list of theoretical masses for matching rows.

  This is the fallback path in _masses_for_label:

  def _masses_for_label(lbl):
      m = comp2masses.get(lbl, [])  # try string lookup first
      if m:
          return m
      return _masses_by_counts(lbl)  # fallback: match by individual counts

  Why the fallback exists: The string lookup (comp2masses) matches canonical label
  strings. If the prediction label and library label canonicalize differently (e.g.,
  one has "K1" and the other has "KDN1"), the string lookup misses. The counts-based
  fallback compares actual numeric values column-by-column, which catches mismatches in
   string representation.

  So: _canon_comp produces a string for comparison. _parse_counts produces numbers for
  comparison. Different strategies for the same goal — finding the theoretical mass for
   a predicted composition.

⏺ The PPM computation loop:                                                            
                                                                                       
  gppm = float(gate_ppm) — Cache the gate threshold as a float outside the loop (avoids
   re-converting every iteration). This is the user-set PPM tolerance (e.g., 10.0).
                                                                                       
  for j, (lbl, pm) in enumerate(zip(labels, obs)): — Iterates through glycan rows in   
  parallel:
  - j — loop index (0, 1, 2, ...), maps to positions in best_ppm and keep arrays       
  - lbl — canonical composition label for this row (e.g., "F1H5N4S2")                  
  - pm — observed precursor mass for this row (e.g., 2076.0234)      
                                                                                       
  Then inside the loop:                                     
                                                                                       
  if not lbl or np.isnan(pm):                               
      continue                          # skip empty labels or missing masses          
                                                                                       
  masses = _masses_for_label(lbl)       # get theoretical masses (string LUT → counts  
  fallback)                                                                            
  if not masses:                                                                       
      continue                          # no library match found                       
   
  m = np.asarray(masses, dtype=float)   # convert to numpy array                       
  diffs = (pm - m) / m * 1e6           # PPM error against each candidate mass
  k = np.nanargmin(np.abs(diffs))       # index of closest match                       
  best_ppm[j] = float(diffs[k])         # store best PPM                               
                                                                                       
  passed = abs(best_ppm[j]) <= gppm     # within tolerance?                            
  keep[j] = passed                      # mark pass/fail                               
                                                                                       
  For each glycan prediction: look up all theoretical masses for that composition,     
  compute PPM error against each, pick the closest one, check if it's within threshold.
   One row at a time, results stored in the pre-allocated arrays.                      

⏺ Correct. The stats breakdown:                                                        
                                                                                       
  - total — all glycan rows (Non-glycan already excluded by the != "Non-glycan" filter)
  - kept — glycans that passed the PPM gate (pred_ok = True)
  - dropped — total - kept                                                             
  - no_lut — dropped because composition not found in library (can't compute PPM at    
  all)                                                                                 
  - no_pm — dropped because observed precursor mass is missing/NaN                     
  - bad_ppm — had both library mass AND observed mass, but PPM exceeded threshold      
                                                                                       
  These three categories (no_lut + no_pm + bad_ppm) should account for all dropped     
  rows. The logger line at 8048-8049 prints all of them in one summary.  
  
⏺ _parse_vec converts probability vectors stored as strings back to Python lists. When 
  proba_vector was stored in pred_df (line ~8124), each cell holds a numpy array. But
  after CSV operations or copying, it might become a string like "[0.05, 0.85, 0.10]". 
  _parse_vec tries JSON parsing, then ast.literal_eval, to recover the list of floats.
                                                                                       
  The if "proba_vector" in pred_df_for_report.columns: block ensures all probability   
  vectors have the same length (padded with zeros if needed) so the reporter can
  process them uniformly. If some rows have 10 classes and others have 8 (shouldn't    
  happen but defensive), they all get padded to target length.

  The else branch handles the case where probabilities were stored as separate columns 
  (proba_F1H5N4S2, proba_Non-glycan, etc.) instead of a single vector column — fills
  missing columns with 0.0.                                                            
                                                            
  This whole block is just data formatting to make summarize_predictions happy — it    
  expects either a proba_vector column with uniform-length lists, or individual proba_*
   columns for every class.   


⏺ The reporter (summarize_predictions) uses it to compute per-prediction confidence    
  metrics — things like:                                                               
                                                                                       
  - max probability — how confident the model was in its top pick                      
  - margin — difference between top-1 and top-2 probabilities (large margin = model is 
  sure)                                                                                
  - entropy — how spread out the probability is (low entropy = confident, high =
  uncertain)                                                                           
                                                            
  These end up in the prediction report as quality indicators per row. A prediction    
  with max_prob=0.95, margin=0.80 is highly confident. One with max_prob=0.35, 
  margin=0.02 is essentially guessing.                                                 
                                                            
  The string conversion is just a storage artifact — the reporter parses it back to    
  numbers via _parse_vec. If proba_vector is already a numpy array (not serialized),
  the harmonize block just pads lengths. Either way, the reporter gets a numeric list  
  to work with.                                             

  For your noted bug about "accepted count is always 0" in the summary — that's likely 
  because the reporter's tau/margin thresholds at line ~8241 (ReportParams(tau=0.60, 
  margin=0.05)) are hardcoded, not reading from your tk variables. Every prediction    
  below 0.60 confidence gets marked as "rejected" in the report regardless of your
  settings.

❯ Oh that's the point I blame GPT (not you) for rejecting everything after it          
  introduced this, and since it's hardcoded I can only get a suspicious report from    
  predcition. So I think final line to edit here is params = ReportParams(tau=0.60,    
  margin=0.05, topk=5, sample_cols=("experiment_title", "sample_name")), correct?      
                                                                                       
⏺ Yes — change to read from tk variables:                                              
                                         
  params = ReportParams(                                                               
      tau=float(thresh_val_var.get()),
      margin=float(margin_val_var.get()),                                              
      topk=5,                                               
      sample_cols=("experiment_title", "sample_name")                                  
  )                                                         
                                                                                       
  Now the reporter uses the same thresholds the user set in Train Settings. If         
  thresholding is disabled (enable_thresh_var = False), you might want to pass
  permissive defaults so nothing gets rejected:                                        
                                                            
  params = ReportParams(
      tau=float(thresh_val_var.get()) if enable_thresh_var.get() else 0.0,
      margin=float(margin_val_var.get()) if enable_thresh_var.get() else 1.0,
      topk=5,                                                                          
      sample_cols=("experiment_title", "sample_name")
  )                                                                                    
                                                            
  With tau=0.0, every prediction is "confident enough." With margin=1.0, no prediction 
  is "indifferent." Result: accepted count matches total count when thresholding is
  off.                                                                                 

⏺ They build the completion messagebox text:                                           
                                                                                       
  - ms1_line — "MS1 feature (expected by model): YES" or "NO" — tells user whether the 
  model was trained with protonatedmass as a feature
  - gate_line — "Precursor gate: ON (10 ppm)" or "Precursor gate: OFF" — tells user    
  whether gating was applied and at what tolerance                                     
  - extra — "\nGated table:\n{path}" if the gated CSV exists, empty string otherwise —
  shows the gated output path only when gate ran and produced results                  
                                                            
  These three lines compose the messagebox.showinfo message so the user sees a summary 
  of what settings were used and where to find the output files. Informational only, no
   logic.           
   


# read_fragment_masses_any
## removed
@Henry: adapted to active extract_fragment_masses
# _ppm_delta
## removed
@Henry: looks like gate helpers leftover, confirmed by Claude Code

# _comp_pat = re.compile(r'(F|H|N|S|G|KDN|A|s|p)\s*([0-9]+)', re.I)
## removed
@Henry: looks like gate helpers leftover, confirmed by Claude Code

# _canonicalize_comp_string
## removed
@Henry: looks like gate helpers leftover, confirmed by Claude Code


# _load_insilico_for_gate
## removed
@Henry: asked Claude Code to confirm if run_prediction already covers this. Seems to be yes.

# _short_id
## line 8253 - 8255
@Henry: probably useful in future, not for now.


# create_unlabeled_from_method
## line 8257 - 8318
@Henry: keep as future CLI integration loader, disabled/dead now


# extract_fragment_masses
## line 8323 - 8363
@Henry: extract features from ion list


# extract_ion_intensities
## line 8366 - 8416 
@Henry: extract/unzip the peaklist/peakintensity to dataframe (big)


# create_unlabeled_dataset
## line 8419 - 8611
@Henry: refactored by Claude Code. Confirmed working on both Legacy exp json and v1 exp json

#  ML Analysis window block
## line 8613 - 8808
@Henry: the status bar is underestimated


#### known issues:
1. Legacy EXP json loaded -> Save will only keep first method exists, others remain Unknown 
and will be discarded if we don't force saving CGA or MAS (the test button that we should remove)
Suggest action: examine the method build, tree behavior (if autosave/update is triggered) on legacy loading.
Need to list relavant functions, probably the context will be clear after we build the user manuals
2. skops error:
loading model trusted must be a list of strings. Before version 0.10 trusted could 
be a boolean, but this is no longer supported, due to a reported CVE-2024-37065.
You can pass the output of get_untrusted_types as trusted to load the data. Be sure
to review the output of the function before passing it as trusted.
After we applied the fix:
  from skops.io import load as sk_load, get_untrusted_types                            
  untrusted = get_untrusted_types(file=model_path)                                     
  model = sk_load(model_path, trusted=untrusted)                                       
Another error:
"persistent IDs in protocol 0 must be
ASCII strings" — which is an old skops file format incompatibility. The joblib
fallback loaded that model fine. New models saved with the current skops version
should load correctly with the get_untrusted_types fix.
Suggest action: fix this before deploying it. load_model_any could be the case

3. minor bug: files["json"] leftover: now should be ["metadata"] or ["method"]. For legacy this field is essential
@Claude Code: files["json"] being a last-writer-wins field.
if you introduce a method resolver that takes a family parameter — something
   like get_active_method(files, family="MAS") that queries _methods by family — then  
  every call site asks for what it needs instead of reading a shared slot. The         
  files["json"] field becomes unnecessary and can be deprecated. v1.1 or 1.2 refactor work.
@Henry: observation in code showing the [_methods] sometimes? pointed back to file["json"].
If we already can switch them to methods (do we have resolver already?) then only legacy calls file["json"],
and we can reassign it to [_methods], safe and clean