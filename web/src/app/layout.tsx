import type { Metadata } from "next";
import "./globals.css";
import ChatWidget from "@/components/ChatWidget";
import { SelectionProvider } from "@/lib/selection-context";
import { ApiKeyProvider } from "@/lib/api-key-context";

export const metadata: Metadata = {
  title: "Steroid Atlas",
  description:
    "Interactive UMAP atlas of steroid-metabolizing enzymes and their small-molecule substrates.",
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en">
      <body className="min-h-screen bg-slate-50 text-slate-900 antialiased">
        <ApiKeyProvider>
          <SelectionProvider>
            {children}
            <ChatWidget />
          </SelectionProvider>
        </ApiKeyProvider>
      </body>
    </html>
  );
}
