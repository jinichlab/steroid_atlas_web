"use client";

import {
  createContext,
  useCallback,
  useContext,
  useState,
  type ReactNode,
} from "react";
import type { ViewKind } from "@/lib/types";

/**
 * Shares the Explore page's current selection (lasso / search / cluster) across
 * things that outlive a single view:
 *  - the globally-mounted ChatWidget (answers "what have I selected");
 *  - the *other* Explore views — switching protein ↔ molecule carries the
 *    selection over and highlights the interacting entities. The bridge is
 *    `proteinAccessions`: for a protein selection it's the selected accessions;
 *    for a molecule selection it's the union of their interacting proteins.
 */

export type SelectionItem = {
  kind: "protein" | "molecule";
  name: string;
  chebi?: string | null;
  smiles?: string | null;
  accession?: string | null;
  gene?: string | null;
  organism?: string | null;
  cluster?: number | string | null;
};

export type SelectionOrigin = "lasso" | "search" | "cluster";

export type SelectionState = {
  view: ViewKind | null;
  origin: SelectionOrigin | null;
  /** Capped list of items (see MAX_SELECTION_ITEMS in the Explore client). */
  items: SelectionItem[];
  /** True count in the selection — may exceed items.length when capped. */
  total: number;
  /** Cross-view bridge (uncapped). Protein selection → the selected accessions;
   *  molecule selection → the union of their interacting protein accessions. */
  proteinAccessions: string[];
  /** Uncapped rowKeys of the selected molecules, for exact same-view restore
   *  (empty for a protein selection). */
  moleculeKeys: string[];
};

export const EMPTY_SELECTION: SelectionState = {
  view: null,
  origin: null,
  items: [],
  total: 0,
  proteinAccessions: [],
  moleculeKeys: [],
};

type SelectionContextValue = {
  selection: SelectionState;
  setSelection: (next: SelectionState) => void;
  clearSelection: () => void;
};

const SelectionContext = createContext<SelectionContextValue | null>(null);

export function SelectionProvider({ children }: { children: ReactNode }) {
  const [selection, setSelectionState] = useState<SelectionState>(EMPTY_SELECTION);

  const setSelection = useCallback(
    (next: SelectionState) => setSelectionState(next),
    [],
  );
  const clearSelection = useCallback(
    () => setSelectionState(EMPTY_SELECTION),
    [],
  );

  return (
    <SelectionContext.Provider
      value={{ selection, setSelection, clearSelection }}
    >
      {children}
    </SelectionContext.Provider>
  );
}

/** Safe outside a provider — returns an empty, inert selection. */
export function useSelection(): SelectionContextValue {
  return (
    useContext(SelectionContext) ?? {
      selection: EMPTY_SELECTION,
      setSelection: () => {},
      clearSelection: () => {},
    }
  );
}
