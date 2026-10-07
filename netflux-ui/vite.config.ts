import react from '@vitejs/plugin-react'
import { defineConfig } from 'vite'

// The FastAPI backend (uvicorn src.api.app:app --port 8000) is proxied under
// /api so the browser never makes a cross-origin request in dev.
export default defineConfig({
  plugins: [react()],
  server: {
    proxy: {
      '/api': { target: process.env.NETFLUX_API ?? 'http://127.0.0.1:8000', changeOrigin: true },
    },
  },
})
