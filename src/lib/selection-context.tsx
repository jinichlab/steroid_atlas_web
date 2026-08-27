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
 * Shares the Explore page's current selection (lasso / search / cluster) with
 * the globally-mounted ChatWidget so the assistant can answer questions about
 * "what I have selected". The Explore client pushes a compact projection of its
 * `pool` into here; ChatWidget reads it and forwards it to /api/chat.
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
};

export const EMPTY_SELECTION: SelectionState = {
  view: null,
  origin: null,
  items: [],
  total: 0,
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
